import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import math





def prior_instant_R_t(t):
    '''
    a = 1
    b = 5
    return np.random.gamma(shape= a, scale = b)
    '''
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




def estimate_R_t(t, pi, incidence_data, w, a = 1, b = 5, n = 1, sample = False):
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
    estimates2 = []
    vals = stats.gamma.ppf([0.025, 0.975], a = alpha, scale = 1/beta)
    #print(vals)

    estimates2.append(alpha/beta)
    if sample:
        for i in range(n):
            estimates.append(np.random.gamma(alpha, 1 / beta))
        return estimates
    return vals, estimates2[0]

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

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


def read_file(file_name):
    """
    :param file_name: reads data from the txt file.
    :return: two arrays denoting days and no_of_cases (I_t)
    """

    txtFile = open(file_name)
    txtFile.readline() #excluding first line
    data = txtFile.readlines()
    days, incidence_data = [],[]
    for line in data:
        line = line.strip()
        line_split = line.split(sep = '\t')
        days.append(int(line_split[0]))
        incidence_data.append(int(line_split[1]))
    return days, incidence_data

def model_epidemic(data_file, prior_rt, mean_si, sd_si, window = 1, uncertain_w = False,  plot_w = False, plot_incidence = False, plot_r_t = False, plot_surface = None, label = ''):
    days, incidence_data = read_file(data_file)
    w = infection_profile(mean_si, sd_si)
    T = days[-1]
    serial_interval_prob = []
    for t in range(T):
        serial_interval_prob.append(w(t))


    if plot_w:

        plot_surface.bar([i + 1 for i in range(T)], serial_interval_prob, width = 0.4, label = label)
        plot_surface.set_xlabel("Days T")
        plot_surface.set_ylabel("Probability of w(t)")



    if plot_incidence:
        print(max(incidence_data))
        plot_surface.plot([i + 1 for i in range(T)], incidence_data, label = label)
        plot_surface.set_xlabel("Days T")
        plot_surface.set_ylabel("Incidence per day")
        plot_surface.grid()


    if plot_r_t:
        n = 1
        for i in range(T):
            if sum(incidence_data[:i + 1]) >= 11 and i > window:
                day = i
                break

        if uncertain_w:
            N = 100
            estimates_of_R_t = []
            for t in range(T):
                posterior_distribution_R_t = np.zeros((N, N))
                for k in range(N):
                    sample_mean = np.random.normal(mean_si, 1.5)
                    sample_SI = sample_mean + 1
                    while sample_SI > sample_mean:
                        sample_SI = np.random.normal(sd_si, 0.5)
                    w = infection_profile(sample_mean, sample_SI)
                    posterior_distribution_R_t[k] = estimate_R_t(t, window, incidence_data = incidence_data, w = w, n = N, sample= True)

                posterior_distribution_R_t = np.ndarray.flatten(posterior_distribution_R_t)
                mean, start, end = mean_confidence_interval(posterior_distribution_R_t)
                median = np.mean(posterior_distribution_R_t)
                estimates_of_R_t.append([(start, end), median])
                print(t, start, end)



            means, start, end = [],[],[]
            for i in range(T):
                m = estimates_of_R_t[i][-1]
                means.append(m)
                start.append(estimates_of_R_t[i][0][0])
                end.append(estimates_of_R_t[i][0][1])
        else:
            estimates_of_R_t = []
            for t in range(T):
                estimates_of_R_t.append(estimate_R_t(t, window, incidence_data=incidence_data, w=w, n=n))

            means, start, end = [], [], []
            for i in range(T):
                m = estimates_of_R_t[i][-1]
                means.append(m)
                start.append(estimates_of_R_t[i][0][0])
                end.append(estimates_of_R_t[i][0][1])




        plot_surface.fill_between([i + 1 for i in range(T)][day:], start[(day):], end[(day):], alpha=0.25)
        #print(estimates_of_R_t.index(max(estimates_of_R_t[day:])))
        plot_surface.plot([i + 1 for i in range(T)][day:], means[day:], label = label)
        plot_surface.set_xlabel("Days T")
        plot_surface.set_ylabel("R_t")


names = ["punjab", "sindh", "GB", "ICT", "KPK", "AJK", "balochistan", "Pakistan", "United Kingdom", "Italy", "Spain"]
names = [i + ".txt" for i in names]

select = [names[-2]]


fig, axs = plt.subplots(1,2)

for i in range(len(select)):
    name = select[i]
    plot_surface = axs[0]
    model_epidemic(name, prior_rt=prior_instant_R_t, window = 7, mean_si= 8.4, sd_si= 3.8, plot_incidence= True, plot_surface= plot_surface)
    plot_surface.axhline(1, ls='--', label = '1')
    plot_surface.legend()
    #plot_surface.set_ylim((0, 5))
    plot_surface.set_title(name[:-4])
    plot_surface = axs[1]
    model_epidemic(name, prior_rt=prior_instant_R_t, window=7, mean_si=8.4, sd_si=3.8, plot_r_t=True,
                   plot_surface=plot_surface)
    plot_surface.axhline(1, ls='--', label='1')
    plot_surface.legend()
    plot_surface.set_ylim((0, 5))
    plot_surface.set_title(name[:-4])

#plt.grid()

"""
fig, axs = plt.subplots(1)
name = "total.txt"
plot_surface = axs

for i in range(len(names[:-1])):
    name = names[i]
    model_epidemic(name, prior_rt=prior_instant_R_t, window = 7, mean_si= 7.5, sd_si= 1, plot_r_t= True, plot_surface= plot_surface, label = name[:-4])
#plot_surface.set_title(name[:-4])
plot_surface.set_ylim((0, 20))

"""
#plt.legend()
plt.tight_layout()
plt.show()