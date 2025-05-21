import numpy as np
import sklearn.cluster as cluster
import sklearn.metrics as metrics
import sklearn.preprocessing
import operator
import scipy.signal as signal

global number_of_points

def normalise(sequence):
	sequence = np.array(sequence).reshape(-1,1)
	scaler = sklearn.preprocessing.MinMaxScaler()
	scaler.fit(sequence)
	norm = scaler.transform(sequence)
	return norm


def calinski_metric(X, labels):
	"""This function is used to calculate Calinski-Harabasz metric"""
	score = metrics.calinski_harabasz_score(X,labels)
	return score


def sil_metric(X, labels):
	"""This function is used to calculate Silhouette metric"""
	score = metrics.silhouette_score(X, labels)
	return score


def davies_metric(X, labels):
	"""This function is used to calculate Davies-Boulding metric"""
	score = metrics.davies_bouldin_score(X, labels)
	return score


def Ncluster_stats(data, start=150, end=500, step=10, plot = False):
#def Ncluster_stats(data, start=1, end=50, step=1, plot = False):
	"""
	This function is used to determine the search range for m-vote
	:param data: pts data
	:param start: search range start
	:param end: search range end
	:param step: search range step
	:return: None
	"""
	"Some empty lists"
	calinski_results = []
	sil_results = []
	davies_results = []

	N = [i for i in range(start, end, step)]

	for n in N:
		#print(n)
		"This part makes clustering for given n"
		estimator = cluster.KMeans(n, n_init=20) #HIGHER n_init = longer computation, better results
		estimator_fit = estimator.fit(pts)
		labels = estimator_fit.labels_
		if len(np.unique(labels)) == len(labels): continue
		
		calinski_results.append(calinski_metric(pts, labels))
		sil_results.append(sil_metric(pts, labels))
		davies_results.append(-davies_metric(pts, labels))

	calinski_results_norm = normalise(calinski_results)
	sil_results_norm = normalise(sil_results)
	davies_results_norm = normalise(davies_results)

	summed = calinski_results_norm + sil_results_norm + davies_results_norm
	detrend = []
	for n in range(1, len(davies_results_norm)):
		detrend.append(summed[n] - summed[n-1])

	highest_cluster = N[np.argmax(np.array(detrend)) + 1]
	del detrend[np.argmax(np.array(detrend))]
	mid_cluster = N[np.argmax(np.array(detrend)) + 1]
	if abs(highest_cluster-mid_cluster) >= number_of_points/20:
		if highest_cluster > mid_cluster:
			search_range = (mid_cluster, highest_cluster)
		else:
			search_range = (highest_cluster, mid_cluster)
	else:
		if highest_cluster > mid_cluster:
			search_range = (mid_cluster - (number_of_points/20), highest_cluster + (number_of_points/20))
		else:
			search_range = (highest_cluster - (number_of_points / 20), mid_cluster + (number_of_points / 20))
	return search_range


def translation(seq):
	"""This function reenumerate sequence from 0 to the end"""
	translation_dict = {}
	cur = 0
	for n in seq:
		if n in translation_dict.keys():
			pass
		else:
			translation_dict[n] = cur
			cur += 1
	for n in range(len(seq)):
		seq[n] = translation_dict[seq[n]]
	return seq


def groups(seq):
	"""This function creates a dictionary with clusters numbers versus points indices
	found in that clusters"""
	community_dict = {}
	seq = list(seq)
	for n in range(len(seq)):
		if seq[n] in community_dict.keys():
			community_dict[seq[n]].append(n)
		else:
			community_dict[seq[n]] = [n]
	return community_dict


def vote(dicts):
	"""This function makes magic
	:param score: weight"""
	import operator
	point_dependencies = {}
	for dic, score in dicts:
		for ind, group in dic.items():
			if ind in point_dependencies.keys():
				pass
			else:
				point_dependencies[ind] = {}
			for dependent in group:
				if dependent in point_dependencies[ind].keys():
					point_dependencies[ind][dependent] += 1 * abs(score)
				else:
					point_dependencies[ind][dependent] = 1 * abs(score)

	majority = {}
	for point in range(number_of_points):
		majority[point] = {}
		for index, cluster in point_dependencies.items():
			if point in cluster.keys():
				majority[point][index] = cluster[point]

	assignment = {}
	for point, cluster in majority.items():
		max_cluster = max(cluster.items(), key=operator.itemgetter(1))[0]
		assignment[point] = max_cluster
	return assignment.values()


def vote_kmeans(pts, max_loops, search_range, last_run = False):
	"""
	KMEANS voting algorithm
	:param pts: data
	:param max_loops: loops limit
	:param search_range: N search range (tuple)
	:return: clustering labels
	"""

	"""Let's initilize some empty variables"""
	labels_stack = []
	calinski_results = []
	sil_results = []
	davies_results = []

	calinski_from_resets = {}

	prev_calin = 0
	prev_sil = 0
	fails = 0
	limiter = 0
	i = 0
	resets = 0

	while i < 20:

		print("i: %s" % i)
		print("Fails: %s" % fails)
		if limiter > max_loops:
			break
		limiter += 1

		#New run
		print("New run")
		ran = np.random.randint(*search_range)
		if last_run:
			estimator = cluster.KMeans(ran, n_init = 100)
		else:
			estimator = cluster.KMeans(ran, n_init = 5)
		estimator_fit = estimator.fit(pts)
		labels = estimator_fit.labels_

		"""Let's calculate metrics for this run"""
		print("Let's calculate metrics for this run")
		calinski_each = calinski_metric(pts, labels)
		sil_each = sil_metric(pts, labels)
		davies_each = davies_metric(pts, labels)

		"""Let's vote!"""
		print("Let's vote!")
		labels = translation(labels) #reenumerate
		labels = groups(labels)
		labels_stack.append((labels, calinski_each)) #second param here is a weight value
		labels = vote(dicts=labels_stack)
		labels = np.array(list(labels)).reshape(-1, 1)

		"""Evaluate current state"""
		print("Evaluate current state")
		calinski_voted = calinski_metric(pts, labels)
		sil_voted = sil_metric(pts, labels)
		davies_voted = davies_metric(pts, labels)

		if prev_calin > calinski_voted or prev_sil > sil_voted:
			print("prev_calin > calinski_voted or prev_sil > sil_voted")
			labels_stack = labels_stack[:-1]
			if last_run:
				print("last_run")
				continue
				i += 1
			else:
				print("not last_run")
				fails += 1
				if fails > 40:
					print("fails > 40")
					resets += 1
					calinski_from_resets[last_succesive_N] = calinski_voted
					calinski_results = []
					if resets >= 15:
						print("resets >= 15")
						max_solution = max(calinski_from_resets.items(), key=operator.itemgetter(1))[0]
						print("Cannot find stable solution, but I reached max Calinski-Harabasz for %s clusters" % str(max_solution))
						print("Fitting %s clusters" % str(max_solution))
						vote_kmeans(pts, max_loops, (max_solution,max_solution+1), True)
						break

					labels_stack = []
					prev_calin = 0
					prev_sil = 0
					i = 0
					fails = 0
					continue
				else:
					print("fails <= 40")
					continue

		last_succesive_N = ran

		i += 1
		fails = 0
		prev_calin = calinski_voted
		prev_sil = sil_voted
		prev_davies = davies_voted

		"""Print some value so we know if clustering have passed and what are the values"""
		print("CALINSKI: %s" % calinski_voted)
		print("SILHOUETTE: %s" % sil_voted)
		print("DAVIES: %s" % davies_voted)


		calinski_results.append(calinski_voted)
		sil_results.append(sil_voted)

		"""STOP if last two clusterings have identical (10e-15) Calinski-Harabasz value as
		two before them"""
		if i >= 4:
			if abs(np.mean(calinski_results[-2:]) - np.mean(calinski_results[-4:-2:])) < 10e-15:
				print("FOUND IT! Had %s clusters" % str(last_succesive_N))
				break

	return labels


#path = """/home/tomowuels/TOMO_OUT/TOMO_2013_05_29_06_14/WORK/mean_points.txt"""
path = """/home/estera/TOMO-LAB_test/EXTERNAL/mvote/mean_points.txt"""

pts = np.loadtxt(path, delimiter=',')
number_of_points = pts.shape[0]
lower_limit = int(round((number_of_points/10)*0.35)*10)
upper_limit = int(round((number_of_points/10)*0.55)*10)

search_range = Ncluster_stats(pts,lower_limit,upper_limit,5)
#search_range = Ncluster_stats(pts,lower_limit,upper_limit,1)
print(search_range)
labels = vote_kmeans(pts, 10000, search_range)

np.savetxt("/home/estera/TOMO-LAB_test/EXTERNAL/mvote/mean_point_mvote_labels", labels, fmt='%.5f', delimiter=' ')
