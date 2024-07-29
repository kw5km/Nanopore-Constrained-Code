# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""

import numpy as np

'''
Franaszek Algorithm for Approx. Eigen Vector 
A <- Graph adjacency matrix
n <- scaling integer
'''
def approx_eigen(A, n):
    for eps in range(1,100):
        y = np.ones((A.shape[0],))*eps
        x = np.zeros_like(y)
        while not np.array_equal(x, y):
            x = y
            Ax = A@x
            nAx = (1/n)*Ax
            y = np.minimum(np.floor(nAx), x)
        if np.max(x)>1: break
    check = np.sum(np.greater_equal(np.matmul(A,x), n*x))
    # print('check:', check)
    # print(f'Stopping on {eps}')
    return x

'''
Reduce Adjecency matrix to sink given approx. eigen vector
A <- Transistion matrix
states <- state names
x <- (A, n)-approx. eigenvector
'''
def sink_graph(A, states, x):

    good_idxs = np.where(x>0)[0]
    new_x = x[good_idxs]
    new_A = A[good_idxs, :][:, good_idxs]
    new_states = states[good_idxs]

    return new_A, new_states, new_x

'''
Splits state i into two states defined by outgoing edge partition E_1, E_2
returns new transition matrix and state names
'''
def add_state(new_A, state, new_states, i, E_1, E_2):
    u_1, c_1 = np.unique(E_1, return_counts=True)
    u_2, c_2 = np.unique(E_2, return_counts=True)   

    # Set state u as u^1
    new_A[i, u_2] = 0    
    new_A[i,u_1] = c_1
    
    # Insert new row and column for u^2 at i+1
    E_2_row = np.zeros_like(new_A[i,:])
    E_2_row[u_2]=c_2
    new_A = np.insert(new_A, i+1, E_2_row , axis=0) 
    
    # Same edges terminate in u^1 and u^2
    E_2_col = new_A[:, i]
    new_A = np.insert(new_A, i+1, E_2_col, axis=1)
    
    # New states names
    new_states = np.insert(new_states, i+1, state)

    return new_A, new_states

'''
iterative Basic x-Consistent Splitting 
A <- Graph adjacency matrix
x <- approx. eigen vector
n <- scaling integer
states <- list of state names corresponding to adjacency matrix row-column indices
'''
def x_split(A, x, n, states, q):
    x_max = np.max(x)
    
    i = 0
    while (x==x_max).any():
        i+=1

        # find all maximal nodes
        ms = np.where(x==x_max)[0]
        # Limit to maximal node meeting Props. 3.3 and 3.4
        vis = [m for m in ms if (x[np.where(A[m, :]>0)]<x_max).any()]

        if len(vis) == 0: break
        vi = vis[0]
         
        state = states[vi]  

        # Edges that origniate from vertex vi
        E_u = np.where(A[vi,:]>0)[0]
        # account for multi-edge state transitions
        E_u = np.hstack([np.tile(E_u[e], int(A[vi, E_u[e]])) for e in range(E_u.size)])

        # Check transitions are correct
        # for e in E_u:
        #     if states[e][:-q] != state[q:]:
                # print(state[q:], states[e][:-q], state, states[e], vi)


        E_u_sorted =E_u# np.concatenate((E_u[v_idx:], E_u[:v_idx])) 

        theta_ms = [np.sum(x[E_u_sorted[:m+1]]) for m in range(E_u.size)]
        p_ms = np.mod(theta_ms, n)
        p_ms = p_ms[:n]

        # Pigeon Hole Case 1
        if 0 in p_ms:
            E_1_idx = np.where(p_ms==0)[0][0]  
            E_1, E_2 = E_u_sorted[:E_1_idx+1], E_u_sorted[E_1_idx+1:] 

        # Pigeon Hole Case 2
        else:
            u, c = np.unique(p_ms, return_counts=True)
            dup = u[c > 1] 
            
            # No zeros no dups? -> WRONG
            if dup.size == 0:
                print('no zero no dupes...','vi:', vi, 'x[i]:', x[i], 'x[E_u]', x[E_u])
                print('p_ms, p_ms.size:', p_ms, p_ms.size, 'dup:', dup)
                break 

            dup_start, dup_end = np.where(p_ms==dup[0])[0][0], np.where(p_ms==dup[0])[0][1] 
            E_1 = E_u_sorted[dup_start+1:dup_end+1]
            E_2 = np.concatenate((E_u_sorted[:dup_start+1], E_u_sorted[dup_end+1:]))

        y_1 = np.sum(x[E_1])*1/n
        y_2 = x[vi] - y_1  
        x[vi] = y_1 
        x = np.insert(x, vi+1, y_2)

        A, states = add_state(A, state, states, vi, E_1, E_2)
    return A, states, x


'''
State-Splitting
A <- Graph adjacency matrix
p <- Cap(S) <= p/q, min out degree: 2^p
q <- power A^q
states <- list of state names corresponding to adjacency matrix row-column indices  
'''
def state_split(A, p, q, states):
    sum_A = np.sum(A, axis=1)
    min_trans = np.min(sum_A[np.nonzero(sum_A)])
    # print('Starting Min Trans:', min_trans)
    
    A_q = np.linalg.matrix_power(A, q)
    x = approx_eigen(A_q, 2**p)
    A_q, states, x = sink_graph(A_q, states, x)
    # print(x.shape)
    # print('x:', x)


    while min_trans<2**p and (np.max(x)>1):
        
        A_q, states, x = x_split(A_q, x,  2**p, states, q)
        
        sum_A_q = np.sum(A_q, axis=1)
        min_trans = np.min(sum_A_q[np.nonzero(sum_A_q)])
        
    #     print('New Min Trans:', min_trans)
    # print('final x:', x)    
    return A_q, states