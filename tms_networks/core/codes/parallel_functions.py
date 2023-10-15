import multiprocessing

def fun(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))

def parallel_map(f, X, nprocs=multiprocessing.cpu_count()):
    q_in = multiprocessing.Queue(1) # '1' is maxsize
    q_out = multiprocessing.Queue() # default maxsize is 2147483647L

    proc = [multiprocessing.Process(target=fun,
                                    args=(f, q_in, q_out)) for _ in range(nprocs)]
    
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    results = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return sorted(results)
