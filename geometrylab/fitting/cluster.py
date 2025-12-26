# -*- coding: utf-8 -*-
from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy.spatial import cKDTree

#import hdbscan

#------------------------------------------------------------------------------

from geometrylab import utilities

#------------------------------------------------------------------------------

__author__ = 'Davide Pellis'

#------------------------------------------------------------------------------

def dbscan(points, radius=None, epsilon=1, min_samples=3):
    '''Density-Based Spatial Clustering of Applications with Noise

    parameters:
        points -- d-dimensional points to cluster, as (n, d) array
    keyword arguments:
        radius -- clustering radius, if None it is computed automatically
        epsilon -- if radius is None, scale the automatic radius
    returns:
        clusters -- (n) array of cluster number starting from 0, -1 is noise
    '''
    N = len(points)
    if radius is None:
        barycenter = np.sum(points, axis=0) / N
        norm = np.linalg.norm(points - barycenter, axis=1)
        radius = (np.max(norm) - np.min(norm))/N**(1/points.shape[1])*epsilon
    tree = cKDTree(points)
    clusters = -np.ones(N, 'i')
    visited = np.full(N, False)
    current_cluster = 0
    visiting = [0]
    while not np.all(visited):
        while len(visiting) > 0:
            next_visiting = []
            for v in visiting:
                if not visited[v]:
                    ball = tree.query_ball_point(points[v], radius)
                    if len(ball) >= min_samples:
                        clusters[v] = current_cluster
                        clusters[ball] = current_cluster
                        a = np.invert(visited[ball])
                        next_visiting.extend(np.array(ball)[a])
                    visited[v] = True
            visiting = np.unique(next_visiting)
        visiting = np.nonzero(np.invert(visited))[0]
        if len(visiting) > 0:
            visiting = [visiting[0]]
        current_cluster += 1
    return clusters

def dbscan_progressive(points, epsilon=.1, size_factor=.1, min_coverage=.7,
                      delta=.1, min_samples=3):
    min_coverage = min(min_coverage, 1)
    min_coverage = max(min_coverage, 0)
    size_factor = min(size_factor, 1)
    size_factor = max(size_factor, 0)
    epsilon = max(.01, epsilon)
    N = len(points)
    if size_factor == 1:
        return np.zeros(N)
    if size_factor == 0:
        return -np.ones(N)
    n_min = max(1, int(N*size_factor))
    min_uncovered = int(max(1, N * (1 - min_coverage)))
    ungrouped = np.full(N, True)
    clusters = -np.ones(N)
    eps = epsilon * 1
    stop = False
    g = 0
    iteration = 0
    while not stop and iteration < 1000:
        C = dbscan(points[ungrouped], epsilon=eps, min_samples=min_samples)
        indices = np.arange(N)[ungrouped]
        for i in range(int(np.max(C)+1)):
            group = np.where(C == i)[0]
            if len(group) > n_min:
                clusters[indices[group]] = g
                g += 1
                ungrouped[indices[group]] = False
        eps *= (1 + delta)
        iteration += 1
        if np.sum(ungrouped) <= min_uncovered:
            stop = True
        if np.sum(ungrouped) < n_min:
            n_min = int(max(1, n_min * (1-delta)))
        if iteration == 1000:
            print('max iteration reached!')
    return clusters

def kmeans(points, number_of_clusters=3, threshold=1e-2,
           return_centroids=False, sample=False):
    Pmax = np.max(points, axis=0)
    Pmin = np.min(points, axis=0)
    if sample:
        pass
    else:
        centers = points[np.random.choice(points.shape[0], number_of_clusters)]
    for i in range(5):
        Ctree = cKDTree(centers)
        D = (Ctree.query(points)[0])**2
        P = D / np.sum(D)
        ind = np.random.choice(points.shape[0], number_of_clusters, p=P)
        centers = points[ind]
    stop = False
    threshold = np.linalg.norm(Pmax-Pmin) * threshold
    while not stop:
        Ctree = cKDTree(centers)
        clusters = Ctree.query(points)[1]
        mask = np.unique(clusters)
        new_centers = np.copy(centers)
        for cluster in mask:
            index = np.where(clusters == cluster)[0]
            C = np.sum(points[index], axis=0) / len(index)
            new_centers[cluster] = C
        if np.max(np.linalg.norm(new_centers - centers, axis=1)) < threshold:
            stop = True
        centers = new_centers
    clusters = utilities.repeated_range(clusters)
    if return_centroids:
        return (clusters, centers)
    else:
        return clusters

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                   Clusterer
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


class Clusterer(object):

    def __init__(self):

        self._points = None

        self._kdtree = None

        self._clusters = None

    @property
    def clusters(self):
        return self._clusters

    def dbscan(self, points, radius=None, epsilon=1, min_samples=3):
        self._points_check(points)
        N = len(points)
        if radius is None:
            barycenter = np.sum(points, axis=0) / N
            norm = np.linalg.norm(points - barycenter, axis=1)
            radius = (np.max(norm) - np.min(norm))/N**(1/points.shape[1])*epsilon
        clusters = -np.ones(N, 'i')
        visited = np.full(N, False)
        current_cluster = 0
        visiting = [0]
        while not np.all(visited):
            while len(visiting) > 0:
                next_visiting = []
                for v in visiting:
                    if not visited[v]:
                        ball = self._kdtree.query_ball_point(points[v], radius)
                        if len(ball) >= min_samples:
                            clusters[v] = current_cluster
                            clusters[ball] = current_cluster
                            a = np.invert(visited[ball])
                            next_visiting.extend(np.array(ball)[a])
                        visited[v] = True
                visiting = np.unique(next_visiting)
            visiting = np.nonzero(np.invert(visited))[0]
            if len(visiting) > 0:
                visiting = [visiting[0]]
            current_cluster += 1
        self._clusters = clusters

    def pdbscan(self, points, epsilon=.1, size_factor=.1, min_coverage=.7,
                delta=.1, min_samples=3):
        min_coverage = min(min_coverage, 1)
        min_coverage = max(min_coverage, 0)
        size_factor = min(size_factor, 1)
        size_factor = max(size_factor, 0)
        epsilon = max(.01, epsilon)
        N = len(points)
        if size_factor == 1:
            return np.zeros(N)
        if size_factor == 0:
            return -np.ones(N)
        n_min = max(1, int(N*size_factor))
        min_uncovered = int(max(1, N * (1 - min_coverage)))
        ungrouped = np.full(N, True)
        clusters = -np.ones(N)
        eps = epsilon * 1
        stop = False
        g = 0
        iteration = 0
        while not stop and iteration < 1000:
            self.dbscan(points, epsilon=eps, min_samples=min_samples)
            C = self._clusters[ungrouped]
            indices = np.arange(N)[ungrouped]
            for i in range(int(np.max(C)+1)):
                group = np.where(C == i)[0]
                if len(group) > n_min:
                    clusters[indices[group]] = g
                    g += 1
                    ungrouped[indices[group]] = False
            eps *= (1 + delta)
            iteration += 1
            if np.sum(ungrouped) <= min_uncovered:
                stop = True
            if np.sum(ungrouped) < n_min:
                n_min = int(max(1, n_min * (1-delta)))
            if iteration == 1000:
                print('max iteration reached!')
        self._clusters = clusters

    def _points_check(self, points):
        if self._points is None:
            changed = True
        elif points.shape[0] != self._points.shape[0]:
            changed = True
        elif points.shape[1] != self._points.shape[1]:
            changed = True
        elif np.linalg.norm(self._points - points) > 1e-10:
            changed = True
        else:
            changed = False
        if changed:
            self._points = points
            self._kdtree = cKDTree(points)



if __name__ == '__main__':
    A = (np.random.random((500,3)) - np.random.random((500,3))) * 2
    #A[:,2] = 0
    #c = kmeans(A, 7)
    c = kmeans(A, 20)
    print(c)

    from geometrylab.vtkplot import Points, view

    pl = Points(A, vertex_data=c, color='Vega20', lut_range='-:+')
    view([pl])