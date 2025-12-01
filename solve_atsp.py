import numpy as np
import scipy.optimize as so
from scipy.optimize import BFGS
import matplotlib.pyplot as plt
import time
import utils
from tqdm import tqdm

import random
import json
import partial_dag_model

def perm_to_tour(perm):
    res = []
    visited = np.zeros([len(perm)])
    c = 0
    while True:
        if visited[c] == 1:
            break
        visited[c] = 1
        res.append(int(c))
        c = perm[c]        
    if np.sum(visited) != len(perm):
        succ = False
    else:
        succ = True

    return succ, res

def assign_to_mat(row, col):
    res = np.zeros([len(row), len(row)])
    res[row, col] = 1.0
    return res


def tour_len(C, tour):
    res = 0.0 
    for i in range(C.shape[0]):
        i_ = tour[i]
        j_ = tour[(i + 1) % len(tour)]
        res += C[i_,j_]
    return res

def solve_TSP(C_, max_loop_cnt=20, visited=None, mu=1.0):
    if visited is None:
        visited = dict()
    nnode = C_.shape[0]
    loop_cnt = 0
    Y = - C_
    X = np.exp(Y)
    
    row, Path = so.linear_sum_assignment(-Y)
    scale_factor = 0.001
    min_scale_factor = 0.0000001
    best_len = 1e20
    best_tour = None 
    while True:
        succ, tour = perm_to_tour(Path)
    
        if succ:
            break
        else:
            start_node = np.random.randint(0, nnode)
            _, tour = partial_dag_model.decodeNNTSP(C_.reshape(-1), Path, start_node)
            nPath = np.zeros(nnode, dtype=np.int32)
            for i_ in range(nnode):
                i = tour[i_]
                j = tour[(i_ + 1) % nnode]
                nPath[i] = j 

            # nPath = partial_dag_model.local_search_2opt(C_.reshape(-1), nPath)
            clength = np.sum(C_[row, nPath])
            if clength < best_len:
                best_len = clength 
                best_tour = nPath

                
        h, gh_ = partial_dag_model.partial_dag_cons(row, Path, 1.0 / 1.05)
        
       
        if json.dumps([int(k) for k in Path]) not in visited.keys():
            visited[json.dumps([int(k) for k in Path])] = False

        max_coeff = np.min((-Y[row, Path]) / (gh_.reshape(nnode, nnode)[row, Path] + 1e-4))
        min_coeff = 1e-5 * max_coeff

        
        default_coeff = scale_factor * max_coeff
        Y = Y - default_coeff * (gh_.reshape(nnode, nnode) + mu * C * X)
        row, diffPath = so.linear_sum_assignment(-Y)
        if np.any(diffPath != Path):
            Path = diffPath
            
            if json.dumps([int(k) for k in diffPath]) in visited.keys():                
                loop_cnt += 1
                if loop_cnt > max_loop_cnt:
                    break
            else:
                loop_cnt = 0
                
            scale_factor *= 0.5
            if scale_factor < min_scale_factor:
                scale_factor = min_scale_factor
        else:
            scale_factor *= 2
        
        X = np.exp(Y - np.max(Y)) 
        X = X / np.sum(X, axis=0)
        X = X / np.sum(X, axis=1)
        Y = 0.99 * Y + 0.01 * np.log(X + 1e-3)

    return succ, row, Path, visited, best_tour

    
nnode = 20
ninstances = 1
all_matrices = np.zeros([ninstances, nnode, nnode])

all_tour_len = np.zeros([ninstances, 4])
rtime = np.zeros([ninstances, 4])

symetric = False 

np.random.seed(1)
for i in tqdm(range(0, ninstances)):
    if symetric:
        D = np.random.uniform(0, 1, [nnode, 2])
        C = np.sqrt(np.sum(np.square(D.reshape(nnode, 2, 1) - D.T.reshape(1, 2, nnode)), axis=1))
    else:
        C = np.random.uniform(0, 1, [nnode, nnode])
    all_matrices[i] = C

np.random.seed(114514)

for i in tqdm(range(0, ninstances)):

    searched = dict()
    C = all_matrices[i]
    nnode = C.shape[0]
    C = C - C * np.eye(nnode)
    C_ = C + np.eye(nnode) * 1e20
    succ = False
    start = time.time()
    succ, row, Path, visited, best_tour = solve_TSP(C_, 20)
    scale = 0
    c_length = 1e20
    if succ:
        Path__ = partial_dag_model.local_search_2opt(C.reshape(-1), Path, symetric)        
        c_length = np.sum(C[row, Path__])
        best_sol = Path
        succ = succ
    else:
        scale = 0.001

    best_tour = partial_dag_model.local_search_2opt(C.reshape(-1), best_tour, symetric)  
    c_length_ = np.sum(C[row, best_tour])
    if c_length_ < c_length:
        c_length = c_length_
        best_sol = best_tour
        

    end = time.time()
    
    print('time = {} length = {}'.format(end - start, c_length_))

