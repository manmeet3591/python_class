from copent import transent
X = np.random.randn(100,1000)
Y = np.random.randn(100,1000)
print(X.shape, Y.shape)

#from itertools import product
def transfer_entropy(x,y):
    return transent(x,y,1)
if __name__ == '__main__':
    iterable = [1, 2, 3, 4, 5]
    iterable = []
    for i in range(10):
        iterable.append((X[i,:], Y[i,:]))
    p = Pool(processes=100)
    data = p.starmap(transfer_entropy, [(X[i,:],Y[i,:]) for i in range(100)])
    #data = p.map(transfer_entropy, iterable)
    p.close()
    print(data)
