# Copyright 2019, Michael Schneider
import logging
import warnings
import numpy as np
from scipy.special import digamma, gammaln
from scipy.special import logsumexp
from sklearn.utils import check_random_state
from sklearn.cluster import KMeans
from sklearn.exceptions import ConvergenceWarning

import sys  # TODO remove
# pylint: disable-msg=C0103
# pylint: disable=too-many-instance-attributes
# pylint: disable-msg=too-many-arguments


class dpgmm():

    """Variational Bayesian estimation of a Gaussian mixture.

    This class allows to infer an approximate posterior distribution over the
    parameters of a Gaussian mixture distribution. The effective number of
    components can be inferred from the data using a truncated stick breaking
    representation. Note that input data is assumed to be 1-dimensional.

    References
    ----------
    1. Blei, David M. and Michael I. Jordan. (2006). "Variational
       inference for Dirichlet process mixtures". Bayesian analysis 1.1
    2. Bishop, C. M. (2006). Pattern recognition and machine learning.
       (New York: Springer).

    Parameters
    ----------
    T : int, defaults to 100.
        An upper bound for the number of mixture components.
    alpha : float | None
        (weight_concentration_prior) The dirichlet concentration of each component on the weight
        distribution (Dirichlet). The higher concentration puts more mass in
        the center and will lead to more components being active, while a lower
        concentration parameter will lead to more mass at the edge of the
        mixture weights simplex. The value of the parameter must be greater
        than 0. If it is None, it's set to ``1. / n_components``.
    means_prior : array-like, shape (self.p_features,) | None
        The prior on the mean distribution (Gaussian).
        If it is None, it's set to 0 for all components.
    means_prior : array-like, shape (self.p_features,) | None
        The prior on the mean distribution (Gaussian).
        If it is None, it's set to 0 for all components.
    a_prior : array-like, shape (self.p_features,) | None
        The prior on the mean distribution (Gaussian).
        If it is None, it's set to 0 for all components.
    b_prior : array-like, shape (self.p_features,) | None
        The prior on the mean distribution (Gaussian).
        If it is None, it's set to 0 for all components.
    tol : float, defaults to 1e-3.
        The convergence threshold. EM iterations will stop when the
        lower bound average gain on the likelihood (of the training data with
        respect to the model) is below this threshold.
    max_iter : int, defaults to 100.
        The number of EM iterations to perform.
    n_init : int, defaults to 1.
        The number of initializations to perform. The result with the highest
        lower bound value on the likelihood is kept. Initialization is
        using KMeans algorithm (random could be chosen instead,
        the code is there)
    init_mode: str, ["kmeans", "random", "both"]
        method to initialize the phi
    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.
    verbose : int, default to 0.
        Enable verbose output. If 0, only warnings are printed. For 1, also
        output from the information level is given. For values greater than 1,
        debugging output is provided.
    verbose_interval : int, default to 10.
        Number of iteration done before the next informational print.

    Attributes
    ----------
    weights_ : array-like, shape (n_components,)
        The weights of each mixture component.
    alpha: float
        (weight_concentration_prior) The higher concentration puts more mass in
        the center and will lead to more components being active, while a lower
        concentration parameter will lead to more mass at the edge of the
        simplex.
    weight_concentration_ : array-like, shape (n_components,)
        The dirichlet concentration of each component on the weight
        distribution (Dirichlet).
    means_ : array-like, shape (n_components, n_features)
        The mean of each mixture component.
    means_prior_ : array-like, shape (n_components, n_features)
        The prior on the mean distribution (Gaussian).
    covariances_ : array-like, shape (n_components, n_features)
        The covariance of each mixture component.
    lbs_ : array-like, shape (n_iter, max_iter)
        Values of lower_bound for all iterations per step
    lower_bound_ : float
        Lower bound value on the likelihood (of the training data with
        respect to the model) of the best fit of inference.
    n_iter_ : int
        Number of step used by the best fit to reach convergence.

    Notes
    --------
    We mostly follow the notation used in (Blei, 2006).
    """
    # TODO finish class description

    def __init__(self, T=120,
                 means_prior=None, means_precision_prior=None,
                 a_prior=None, b_prior=None,
                 tol=1e-3, max_iter=500, n_init=3, init_method="both",
                 random_state=None, verbose=1, verbose_interval=50):

        self.T = T  # Truncation parameter of infinite sum

        # configuration parameters
        self.tol = tol
        self.init_method = init_method
        self.max_iter = max_iter
        self.n_init = n_init
        self.lbs_ = np.zeros((self.n_init, self.max_iter))
        self.lbs_.fill(np.nan)

        # prior values
        self.means_prior = means_prior # m_k
        self.means_precision_prior = means_precision_prior  # \beta_k

        self.a_prior = a_prior
        self.b_prior = b_prior

        self.s_1_prior = 1.
        self.s_2_prior = 1.

        # model output variables
        #  self.means_precision_ # \beta_k
        #  self.covariance_ = None  # self.a_ / self.b_
        #  self.weight_concentration_ #

        self.random_state = random_state
        self.verbose = verbose
        self.verbose_interval = verbose_interval

        # state variables
        self.isFitted_ = False
        self.lower_bound_ = None

        # Setup logging
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s-%(levelname)s-%(message)s',
                            datefmt='%m/%d/%Y_%H:%M:%S')

        self.logger = logging.getLogger('mixture-model')
        if self.verbose <= 0:
            self.logger.setLevel(logging.WARNING)
        elif self.verbose == 1:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.DEBUG)

    #################### UPDATE EQUATIONS ####################

    def compute_eta_z(self):
        """"
        Returns
        -------
        eta_z : float

        see Blei, p.129 and supplementary materials
        E_{q} [\log V_t] + \sum_{i=1}^{t-1} E_{q} [\log (1-V_i)]...
        """
        digamma_sum = digamma(self.weight_concentration_[0] +
                              self.weight_concentration_[1])
        digamma_a = digamma(self.weight_concentration_[0])
        digamma_b = digamma(self.weight_concentration_[1])

        eta_z = digamma_a - digamma_sum + \
            np.hstack((0, np.cumsum((digamma_b - digamma_sum)[:-1])))

        return eta_z

    def compute_eta_x(self, X):
        """
        Parameters
        ----------
        X : array-like, shape (n_samples, p_features)

        Returns
        -------
        eta_x: float

        see supplementary materials
        """
        eta_x = -.5 * (np.log(2*np.pi) - digamma(self.a_)
                    + np.log(self.b_)
                    + (1./self.means_precision_)
                    + (self.a_/self.b_) * np.square(self.means_ - X))

        return eta_x

    def e_step(self, X):
        """ update phi / equivalent to E step.

        Parameters
        ----------
        X : array-like, shape (n_samples, p_features)

        Returns
        -------
        log_resp: array, shape (n_samples, T)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.
        """
        eta_z = self.compute_eta_z()
        eta_x = self.compute_eta_x(X)

        phi_nk = (eta_z + eta_x - 1)
        assert(phi_nk.shape == (X.shape[0], self.T))

        # normalize phi_nk to probability
        log_sum = logsumexp(phi_nk, axis=1)
        with np.errstate(under='ignore'):
            # ignore underflow
            log_resp = phi_nk - log_sum[:, np.newaxis]

        return log_resp

    def update_alpha(self):

        dg_sum = digamma(self.weight_concentration_[0]
                         + self.weight_concentration_[1])
        dg0 = digamma(self.weight_concentration_[0]) - dg_sum
        dg1 = digamma(self.weight_concentration_[1]) - dg_sum

        self.w_1 = self.s_1_prior + self.T - 1
        self.w_2 = self.s_2_prior - np.sum((dg1 - dg_sum)[:-1])


    def update_V(self, resp):
        """ update V.

        Parameters
        ----------
        resp : array-like, shape (n_samples, p_features)
            posterior probabilities (or responsibilities) of
            the point of each sample in X.

        """
        nk = np.sum(resp, axis=0)

        gamma_t1 = 1 + nk
        gamma_t2 = (self.w_1 / self.w_2) + \
            np.hstack((np.cumsum(nk[::-1])[-2::-1], 0))
        #  gamma_t2 = self.alpha + \
            #  np.hstack((np.cumsum(nk[::-1])[-2::-1], 0))

        # Dirichlet process weight_concentration is a tuple
        # containing the two parameters of the beta distribution
        # see Blei, (Eq. 18, 19) on p. 129
        self.weight_concentration_ = (gamma_t1, gamma_t2)

    def update_mu(self, X, resp):
        """ update \mu.

        Parameters
        ----------
        X : array-like, shape (n_samples, p_features)

        resp : array-like, shape (n_samples, p_features)
            posterior probabilities (or responsibilities) of
            the point of each sample in X.

        """
        nk = np.sum(resp, axis=0)
        self.means_precision_ = self.means_precision_prior + nk

        self.means_ = self.means_prior[:, 0] * self.means_precision_prior + np.dot(resp.T, X)[: , 0]
        self.means_ /= self.means_precision_prior + nk

    def update_Lambda(self, X, resp):
        """ update \Lambda.

        Parameters
        ----------
        X : array-like, shape (n_samples, p_features)
        resp : array-like, shape (n_samples, p_features)
            posterior probabilities (or responsibilities) of
            the point of each sample in X.

        """
        #  print(resp.shape)
        #  print(X.shape)
        #  print(self.means_.shape)
        #  sys.exit(2)
        self.a_ = 1 + .5 * np.sum(resp, axis=0)
        self.b_ = 1 + .5 * np.sum(resp * np.square(self.means_ - X), axis=0) + \
            .5 * self.means_precision_prior * np.square(self.means_ - self.means_prior[:, 0])

    def m_step(self, X, log_resp):
        """M step

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        log_resp : array-like, shape (n_samples, n_components)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.
        """
        # TODO more numerically stable?
        resp = np.exp(log_resp)
        #  resp = resp + np.finfo(resp.dtype).eps
        #  resp = resp / np.sum(resp, axis=1)[:, np.newaxis]

        self.update_V(resp)
        self.update_alpha()
        self.update_mu(X, resp)
        self.update_Lambda(X, resp)

        self.covariance_ = self.b_ / self.a_
        if np.any(np.less_equal(self.covariance_, 0.0)):
            raise ValueError("Fitting the mixture model failed "
                             "because some components have ill-defined "
                             "empirical covariance (for instance caused "
                             "by singleton or collapsed samples). Try to "
                             "decrease the number of components, or increase "
                             "reg_covar.")

    #################### END UPDATE EQUATIONS ####################

    #################### VARIATIONAL LOWER BOUND ####################
    def compute_variational_bound(self, X, log_resp):
        """Estimate the lower bound of the model.

        The lower bound on the likelihood (of the training data with respect to
        the model) is used to detect the convergence and has to decrease at
        each iteration.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        log_resp : array, shape (n_samples, T)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.

        Returns
        -------
        lower_bound : float
            value of the variational lower bound
        """
        lower_bound = 0

        # TODO inluce term for lower bound (alpha)

        resp = np.exp(log_resp)
        #  resp = resp + np.finfo(np.float32).eps # probably necessary and fine TODO
        dg_sum = digamma(self.weight_concentration_[0]
                         + self.weight_concentration_[1])
        dg0 = digamma(self.weight_concentration_[0]) - dg_sum
        dg1 = digamma(self.weight_concentration_[1]) - dg_sum

        # \alpha ##
        # E_q[log p(V | 1, alpha)]
        lp_alpha = -1 * (self.w_1 / self.w_2)
        lq_alpha = -self.w_1 + np.log(self.w_2) - gammaln(self.w_1) - (1 - self.w_1) * digamma(self.w_1)
        lower_bound += lp_alpha - lq_alpha

        # V ##
        # E_q[log p(V | 1, alpha)]
        alpha = self.w_1 / self.w_2
        lp_V = self.T * (gammaln(1 + alpha) - gammaln(alpha)) \
            + (alpha - 1) * np.sum(dg1)
        # E_q[log q(V | gamma1, gamma2)]
        lq_V = np.sum(gammaln(self.weight_concentration_[0] +
                              self.weight_concentration_[1])
                      - gammaln(self.weight_concentration_[0])
                      - gammaln(self.weight_concentration_[1])
                      + (self.weight_concentration_[0] - 1) * dg0
                      + (self.weight_concentration_[1] - 1) * dg1)
        self.logger.debug("LB(V     ): %.5f (%.5f, %.5f)", lp_V - lq_V, lp_V, lq_V)
        lower_bound += (lp_V - lq_V)

        # \mu, \Lambda ##
        lpq_ml = 0
        lpq_ml += .5 * np.sum(np.log(self.means_precision_prior / self.means_precision_) -
                        (self.means_precision_prior / self.means_precision_) + 1
                        - self.means_precision_prior * (self.a_ / self.b_) * np.square(self.means_ -
                                                                                       self.means_prior[:, 0]))
        lpq_ml += np.sum(self.a_ - np.log(self.b_) + gammaln(self.a_) + (1 - self.a_) * digamma(self.a_))
        lpq_ml += np.sum(self.a_prior * np.log(self.b_prior) + (self.a_prior - 1) * (digamma(self.a_) - np.log(self.b_))
                         - self.b_prior * (self.a_ / self.b_) - gammaln(self.a_prior))

        self.logger.debug("LB(normal): %.5f", lpq_ml)
        lower_bound += lpq_ml

        # Z ##
        # Blei, p.129, E_q [\log (z_n | V)]
        dg_cumsum = np.hstack((0, np.cumsum(dg1[:-1])))

        lp_Z = 0
        # E_q[log p(Z | V)]
        for t in range(self.T):
            # TODO debug should we include T in first sum
            lp_Z += np.sum(resp[:, t] * dg_cumsum[t], axis=0)
            lp_Z += np.sum(resp[:, t] * dg0[t], axis=0)

        # E_q[log q(Z)]
        # neg. Entropy: \sum_{i} \phi_{n,k} \log(\phi_{n,k}
        # TODO nan to num necessary
        lq_Z = np.sum(np.nan_to_num(log_resp * resp))

        self.logger.debug("LB(Z     ): %.5f (%.5f, %.5f)", lp_Z - lq_Z, lp_Z, lq_Z)
        lower_bound += (lp_Z - lq_Z)

        # TODO check -> dot product
        lp_X = -.5 * np.sum(np.sum(resp * (np.log(2 * np.pi) -
                                   digamma(self.a_) + np.log(self.b_) +
                                   (1./self.means_precision_) +
                                   ((self.a_ / self.b_) *
                                   np.square(self.means_ - X))),
                                   axis=1), axis=0)
        self.logger.debug("LB(X  ): %.5f", lp_X)
        lower_bound += lp_X

        return lower_bound

    #################### END VARIATIONAL LOWER BOUND ####################

    def _get_parameters(self):
        return (self.weight_concentration_,
                self.means_,
                self.covariance_)

    def _set_parameters(self, params):
        (self.weight_concentration_, self.means_, self.covariance_) = params

        # Weights computation
        weight_dirichlet_sum = (self.weight_concentration_[0] +
                                self.weight_concentration_[1])
        tmp = self.weight_concentration_[1] / weight_dirichlet_sum
        self.weights_ = (
            self.weight_concentration_[0] / weight_dirichlet_sum *
            np.hstack((1, np.cumprod(tmp[:-1]))))
        self.weights_ /= np.sum(self.weights_)

    def initialize_parameters(self, X, random_state):
        #  batch_size, n_samples, p_features = X.shape
        # TODO
        n_samples, p_features = X.shape

        # Make sure X is in the right format
        assert(n_samples > p_features)

        self.p_features = p_features  # dimension of problem
        assert(self.p_features == 1)

        # Choose default priors
        #  if self.alpha is None:
            #  self.logger.info("alpha (weight concentration prior) set automatically")
            #  self.alpha = 1./self.T

        if self.means_prior is None:
            self.logger.info("Mean prior set automatically")
            self.means_prior = np.zeros((self.T,))

        if self.means_precision_prior is None:
            self.logger.info("Mean precision prior set automatically")
            self.means_precision_prior = np.ones((self.T,))

        if self.a_prior is None:
            self.logger.info("Covariance gamma shape prior set automatically")
            self.a_prior = np.ones((self.T,))

        if self.b_prior is None:
            self.logger.info("Covariance gamma rate prior set automatically")
            self.b_prior = np.ones((self.T,))

        # TODO check KMEANs, there might be a bug in there
        # resp initialisation
        if (self.init_method == "both" and self.iter < 2) or self.init_method.lower() == "kmeans":
            #  Kmeans initialisation
            self.logger.debug("Kmeans initialization at step ", self.iter)
            rand_init = random_state.randint(2**32-1, size=self.n_init)
            resp = np.zeros((n_samples, self.T))
            label = KMeans(n_clusters=self.T, n_init=1,
                           random_state=rand_init[self.iter]).fit(X).labels_
            resp[np.arange(n_samples), label] = 1
            """
               alternative is random initialisation
               but requires far more iterations!
            """
        else:
            self.logger.debug("Random initialization at step %d", self.iter)
            resp = random_state.rand(n_samples, self.T)
            resp = resp + np.finfo(np.float32).eps # probably necessary and fine TODO
            resp = np.divide(resp, np.sum(resp, axis=1)[:, np.newaxis])

        self.w_1 = 1
        self.w_2 = 1
        # TODO test some stuff to improve numerical stability
        log_resp = np.log(resp) #+ np.finfo(np.float32).eps)
        self.m_step(X, log_resp)
        self.logger.debug("Model initialized - step %d", self.iter)

    def fit(self, X):
        # TODO include option to do this with or without ELBO
        # get method to compute ELBO for final model
        max_lower_bound = -np.infty
        random_state = check_random_state(self.random_state)

        for i_init in range(self.n_init):

            self.logger.info("Iteration: %d", i_init)

            isConverged = False
            self.iter = i_init
            self.initialize_parameters(X, random_state)
            lower_bound = 0
            prev_lower_bound = -np.infty

            for i_iter in range(self.max_iter):
                self.logger.debug("Iteration:%03d - Step:%04d - "
                                 "LB:%.5f", i_init, i_iter,
                                 lower_bound)

                log_resp = self.e_step(X)
                a = np.exp(log_resp)

                #sys.exit(2)
                self.m_step(X, log_resp)
                lower_bound = self.compute_variational_bound(X, log_resp)

                self.lbs_[i_init, i_iter] = lower_bound
                delta = lower_bound - prev_lower_bound

                if lower_bound <= prev_lower_bound:
                    warnings.warn("Lower bound error")
                    self.logger.warn("Issue with convergence: "
                                     "step: %d, "
                                     "Lower bound: %.5f, "
                                     "Prev lower bound: %.5f, "
                                     "Delta: %.5f", i_iter, lower_bound,
                                     prev_lower_bound, delta)

                diff = lower_bound - prev_lower_bound
                # This has to hold true for all diff,
                # bound is #TODO
                #  assert(diff > 0)
                prev_lower_bound = lower_bound

                if abs(delta) < self.tol:
                    isConverged = True
                    self.logger.info("Converged with delta %.5f", delta)
                    break

                if i_iter % self.verbose_interval == 0:
                    self.logger.info("Iteration:%03d - Step:%04d - "
                                     "LB:%.5f", i_init, i_iter,
                                     lower_bound)

            if lower_bound > max_lower_bound:
                self.logger.debug("Update solution in step %d", self.iter)
                max_lower_bound = lower_bound
                best_params = self._get_parameters()
                best_n_iter = i_iter
                self.lower_bound_ = lower_bound
                self.isFitted_ = True

        if not isConverged:
            warnings.warn('Initialization %d did not converge. '
                          'Try different init parameters, '
                          'or increase max_iter, tol '
                          'or check for degenerate data.'
                          % (i_init + 1), ConvergenceWarning)

        self._set_parameters(best_params)
        self.n_iter_ = best_n_iter

        return (best_params, best_n_iter)

    # TODO rewrite and refactor below
    def score_samples(self, X):
        """Compute the weighted log probabilities for each sample.

        Parameters
        ----------
        X : array-like, shape (n_samples, self.p_features)
            List of self.p_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        log_prob : array, shape (n_samples,)
            Log probabilities of each data point in X.
        """
        assert(self.isFitted_)
        _, _, weighted_log_prob = self.e_step(X)
        return logsumexp(weighted_log_prob, axis=1)

    def score(self, X):
        """Compute the per-sample average log-likelihood of the given data X.
        Parameters
        ----------
        X : array-like, shape (n_samples, n_dimensions)
            List of self.p_features-dimensional data points. Each row
            corresponds to a single data point.
        Returns
        -------
        log_likelihood : float
            Log likelihood of the Gaussian mixture given X.
        """
        assert(self.isFitted_)
        return self.score_samples(X).mean()

    def predict(self, X):
        """Predict the labels for the data samples in X using trained model.
        Parameters
        ----------
        X : array-like, shape (n_samples, self.p_features)
            List of self.p_features-dimensional data points. Each row
            corresponds to a single data point.
        Returns
        -------
        labels : array, shape (n_samples,)
            Component labels.
        """
        assert(self.isFitted_)
        _, _, weighted_log_prob = self.e_step(X)
        return weighted_log_prob.argmax(axis=1)

    def predict_proba(self, X):
        """Predict posterior probability of each component given the data.

        Parameters
        ----------
        X : array-like, shape (n_samples, self.p_features)
            List of self.p_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -, betaln------
        resp : array, shape (n_samples, n_components)
            Returns the probability each Gaussian (state) in
            the model given each sample.
        """
        assert(self.isFitted_)
        _, log_resp, _ = self.e_step(X)
        return np.exp(log_resp)

    def sample(self, n_samples=1):
        """Generate random samples from the fitted Gaussian distribution.

        Parameters
        ----------
        n_samples : int, optional
            Number of samples to generate. Defaults to 1.

        Returns
        -------
        X : array, shape (n_samples, self.p_features)
            Randomly generated sample
        y : array, shape (nsamples,)
            Component labels
        """
        assert(n_samples >= 1)
        assert(self.isFitted_)

        rng = check_random_state(self.random_state)
        n_samples_comp = rng.multinomial(n_samples, self.weights_)

        X = np.vstack([mean + rng.randn(sample, self.p_features) * np.sqrt(covariance)
                       for (mean, covariance, sample) in
                       zip(self.means_, self.covariance_, n_samples_comp)])

        y = np.concatenate([j * np.ones(sample, dtype=int)
                            for j, sample in enumerate(n_samples_comp)])

        return (X, y)


def fit_mixture(X, T,
                random_state=2018, verbose=1,
                n_init=3, max_iter=1000,
                tol=1e-3):

    np.random.seed(random_state)

    alpha = 0.1
    means_prior = create_prior(T, 1, offset=0.0, scale=0.0)
    means_precision_prior = create_prior(T, 1, offset=.1, scale=0)[:, 0]
    a_prior = np.ones((T, ))
    b_prior = np.ones((T, ))

    model = dpgmm(T=T, alpha=alpha,
                  means_prior=means_prior,
                  means_precision_prior=means_precision_prior,
                  a_prior=a_prior, b_prior=b_prior,
                  random_state=random_state, verbose=verbose,
                  n_init=n_init, max_iter=max_iter,
                  tol=tol)

    model.fit(X)

    # collect results
    #  weights = model.weights_
    #  xi = model.xi_
    #  peaks = model.means_ / xi

    # get maximum index, such that cumulative sum of weights is <= .99
    #  idx = np.argsort(weights)[::-1]
    #  max_V = np.sum(np.cumsum(weights[idx]) <= 0.999)
    #  relevant_peaks = np.sort(peaks[idx][:max_V])

    # also return cluster means for visualization purposes
    import pandas as pd
    d = {'means': model.means_,
         'covariance': model.covariance_,
         'weights': model.weights_,
         'index': np.array([i for i in range(model.means_.shape[0])])
        }

    return pd.DataFrame(data=d)
