import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import math





def prior_instant_R_t(t):
    return  2.5
    if t <= 15:
        return 2.5
    return 0.7

def estimate_I_t(t, instant_R_t, incidence_data, w):
    '''
    :param t: t (in days) at which to estimate I_t
    :param instant_R_t: a function of the form f(t) which returns the reproduction number at a given time t.
    :param incidence_data: a list of integers. indcidence_data[s] represents the number of incident cases on day s.
    :param w: a probability density function describing the infectious profile with signature f(t) where is an integer. w(s) represents the probability of spreading the infection on day s.
    :return: prediction of the number of incidence cases on day t.

    Follwing appendix 1 of cori 2013. I_t follows a poisson distribution. (See eq. 1 of equations page)
    '''

    #summation = 0
    #for s in range(t):
    #    summation += (incidence_data[s] * w(s))

    summation = lambda_t(t, incidence_data, w)


    mean = prior_instant_R_t(t) * summation #mean of poisson distribution

    estimate = np.random.poisson(mean) #sample
    return estimate

def lambda_t(t, incidence_data, w):
    """
    :param indidence_data: a list of integers. indcidence_data[s] represents the number of incident cases on day s.
    :param w: a probability density function describing the infectious profile with signature f(t) where is an integer. w(s) represents the probability of spreading the infection on day s.
    :return: a float representing the summation

    Following appendex 1 of cori 2013. See (eq.3 of equations)
    """
    summation = 0
    for s in range(1, t+1):
        summation += (incidence_data[t-s] * w(s))

    return summation




def estimate_R_t(t, pi, incidence_data, w, a = 1, b = 5, n = 1):
    """
    :param t: t (in days) at which to re_estimate R_t
    :param pi: the size of the window
    :param incidence_data: a list of integers. indcidence_data[s] represents the number of incident cases on day s.
    :param w: a probability density function describing the infectious profile with signature f(t) where is an integer. w(s) represents the probability of spreading the infection on day s.
    :param a, b: constants to be used in gamma distribution of a.
    :return: an estimate of Reproduction number on day t.
    Follwing appendix 1 of cori 2013. I_t follows a poisson distribution. (See eq. 2 of equations page)
    """

    summation = 0
    summation_lambdas = 0
    for s in range(t - pi, t):
        summation += (incidence_data[s])
        summation_lambdas += lambda_t(s, incidence_data, w)

    alpha = a + summation
    beta = (1/b) + summation_lambdas

    estimates = []
    for i in range(n):
        estimates.append(np.random.gamma(alpha, 1/beta))
    estimate2 = alpha/beta
    return estimates

def infection_profile(mean, std_deviation):
    """
    :param mean: mean of the gamma distribution that the infection profile follows.
    :param std_deviation: std_deviation of the gamma distribution that the infection profile follows.
    :param s: day at which to determine probability.
    :return: descritized probability of the infection profile at day s.

    Following shifted gamma distribution on appendix 11
    """

    beta = mean / (std_deviation ** 2) #since, mean = alpha/beta and variance = alpha/(beta^2)
    alpha = mean * beta

    probs = {}
    def prob(s):
        int_s = int(s)
        assert np.allclose([int_s], [s])  # making sure s is an integer
        s += 1
        if s in probs:
            return probs[s]
        else:
            prob = (np.convolve(s,stats.gamma.cdf(s, a = alpha, scale = 1/beta))) + (np.convolve((s - 2), stats.gamma.cdf(s - 2, a = alpha, scale = 1/beta)) ) - \
                   (np.convolve(2,  np.convolve((s - 1) , stats.gamma.cdf(s - 1, a = alpha, scale = 1/beta)))) + \
                   ((alpha * 1/beta)*(2*(stats.gamma.cdf(s-1, a = alpha + 1, scale = 1/beta))- stats.gamma.cdf(s - 2, a = alpha + 1, scale = 1/beta) - \
                                      stats.gamma.cdf(s, a = alpha + 1, scale = 1/beta)))
        probs[s] = prob[0]
        return prob[0]
    return prob

if __name__ == "__main__":
    """Following appendix 6"""
    T = 50
    pandemics = 1
    window_size = 1
    N = 100

    incidence_across_simulations = [[] for i in range(1, T+1)]
    reproduction_across_simulations = [[] for i in range(window_size, T+1)]
    min_day = []

    means = []
    start = []
    end = []
    days = []



    for i in range(pandemics):
        mean = 8.4
        std_deviation = 3.8
        incidence_data = [0,10]
        w = infection_profile(mean, std_deviation)
        reproduction_number = [0]
        day_found = False
        probs = []

        for t in range(10):
            probs.append(w(t))

        for t in range(2, T+1):

            i_t_estimates = []
            r_t_estimates = np.zeros((N, N))
            for iter_num in range(N):
                print(t,iter_num)
                r_t_estimates_mean_constant = []
                current_incidence = incidence_data.copy()
                mean_si = np.random.normal(loc = mean, scale = 1.5)
                std_si = math.inf

                while std_si > mean_si:
                    std_si = np.random.normal(loc = std_deviation, scale = 0.5)
                w = infection_profile(mean_si, std_si)

                current_incidence.append(estimate_I_t(t, instant_R_t= prior_instant_R_t, w = w, incidence_data=current_incidence))

                if t > window_size:
                    r_t_estimates[iter_num] = (estimate_R_t(t, window_size, incidence_data=current_incidence, w=w, n = N))
                if not day_found and np.sum(current_incidence) > 12:
                    min_day.append(t)
                    day_found = not day_found

                i_t_estimates.append(current_incidence[-1])

            incidence_data.append(np.mean(i_t_estimates))

            confidence = 0.95
            data = np.ndarray.flatten(r_t_estimates)
            m = np.mean(data)
            std_err = stats.sem(data)
            h = std_err * stats.t.ppf((1 + confidence) / 2, (N*N) - 1)
            start.append(m - h)
            end.append(m + h)
            days.append(t)
            if t > 7:
                plt.scatter([t] * (N * N), np.ndarray.flatten(r_t_estimates), color = "black", marker='x', s = 0.1)

    print(days)
    print(start)
    plt.fill_between(days[max(min_day):], start[max(min_day):], end[max(min_day):], color="blue", alpha=1)



    plt.xlabel("Time (Days)")
    plt.ylabel("Instant R_t")
    plt.show()





    """
        for i in range(len(incidence_data[1:])):
            incidence_across_simulations[i].append(incidence_data[i+1])

        for i in range(len(reproduction_number)):
            reproduction_across_simulations[i].append(reproduction_number[i])
        """


    """ 
        plt.scatter(range(10), probs)
        plt.xlabel("Serial Interval (Days)")
        plt.ylabel("Probability")
        plt.show()
        """
        #print(i)
        #plt.scatter(range(10, T), reproduction_number[10:], color = "black", marker='x')

        #print(min(incidence_data), max(incidence_data), np.mean(incidence_data))

    #means = [np.mean(i) for i in incidence_across_simulations]
    #plt.plot(range(T), means,  color = "red")

    means = []
    start = []
    end = []
    days = []
    #print(min_day)
    #print(max(min_day))
    """
    for i in range(len(reproduction_across_simulations)):
        confidence = 0.95
        data = reproduction_across_simulations[i]
        n = len(data)
        m = np.mean(data)
        means.append(m)
        days.append(i)
        std_err = stats.sem(data)
        h = std_err * stats.t.ppf((1 + confidence) / 2, n - 1)
        start.append(m - h)
        end.append(m + h)
    """



    """
    each_reproduction = np.transpose(reproduction_across_simulations)
    for i in range(len(each_reproduction)):
        plt.plot(days[max(min_day):], each_reproduction[i][max(min_day):], color = "black", linewidth = 0.25)

    plt.fill_between(days[max(min_day):], start[max(min_day):], end[max(min_day):], color="blue", alpha = 1)
    #plt.plot(days[max(min_day):], means[max(min_day):], color = "black")

    plt.xlabel("Time (Days)")
    plt.ylabel("Instant R_t")
    plt.show()

    """











