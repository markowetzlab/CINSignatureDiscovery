# Copyright, 2019, Michael Schneider
import matplotlib.pyplot as plt
import numpy as np

def plot_solution(X, Y):
    # the histogram of the data
    bins = np.linspace(0, np.max(X) + 10, 1000)
    n, bins, patches = plt.hist(X, bins, density=True, facecolor='g', alpha=0.75)

    n, bins, patches = plt.hist(Y, bins, density=True, facecolor='r', alpha=0.25)

    #  alpha, beta = weights
    #  modes = (alpha - 1) / (alpha + beta - 2)
    #  print(modes)
#
    #  for xc in modes:
        #  plt.axvline(x=xc)

    #  plt.xlabel('Smarts')
    #  plt.ylabel('Probability')
    #  plt.title('Histogram of IQ')
    #  plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #  plt.axis([40, 160, 0, 0.03])
    #  plt.grid(True)
    plt.show()

def create_prior(K, D, offset, scale, batch=None):
    prior = offset + (scale * np.reshape(np.repeat(np.array(range(K)), D, axis=0), (-1, D)))
    if batch:
        prior = np.repeat(prior[np.newaxis, :, :], batch, axis=0)

    return prior


# Generate simulated data
def build_dataset(N, K, D, prob, zero_prob, dist):
    """Build a mixture dataset for testing purposes with
        an arbitrary generator function distExample function with types documented in the docstring.
        TODO finish

       Args:
            N (int): number of elements
            K (int): number of mixture components
            D (int): dimensionality of each component
            prob (array-like), shape

       Returns:
            X: The return value. True for success, False otherwise.
            label:
    """
    label = np.random.choice(range(K), size=N, replace=True, p=prob)
    assert(prob.shape[0] == K)

    X = np.zeros((N, D))
    for i in range(N):
        new_element = dist()
        assert(new_element.shape == (D, K))
        X[i, :] = new_element[:, label[i]]

    # dropout events
    idx_zero = np.random.choice(range(N), size=int(zero_prob*N), replace=False)
    for i in idx_zero:
        X[i, :] = np.zeros((D, 1))

    return X, label

