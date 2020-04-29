import numpy as np
from scipy import stats, spatial
import matplotlib.pyplot as plt
import math
import datetime
import matplotlib.dates as mdates






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

def predict_I_t(t, r_t, incidence_data, w):
    summation = lambda_t(t, incidence_data, w)

    prediction = r_t * summation
    conf = stats.poisson.interval(0.95, prediction)
    return conf, prediction

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


def truncated_normal(mean, sd, lower, upper, n = 1):
    a, b = (lower - mean) / sd, (upper - mean) / sd
    return stats.truncnorm.rvs(a,b, loc = mean, scale = sd, size = n)


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

def plot_serial_interval(t, probability, label = ''):
    t = list(t)
    plot_surface.bar(t, probability, width=0.4, label=label)
    plot_surface.set_xlabel("Days T")
    plot_surface.set_ylabel("Probability of w(t)")


def plot_incidence_data(t, data, label = ''):
    plot_surface.bar(list(t), data, label=label)
    plot_surface.set_xlabel("Days T")
    plot_surface.set_ylabel("Incidence per day")
    plot_surface.grid()

def plot_estimates_r_t(t, starts, ends, means, label = ''):
    assert len(t) == len(starts) == len(ends) == len(means)
    plot_surface.fill_between(list(t), starts, ends, alpha = 0.1, color = "black")
    # print(estimates_of_R_t.index(max(estimates_of_R_t[day:])))
    plot_surface.plot(list(t), means, label=label, c='k')
    plot_surface.set_xlabel("Days T")
    plot_surface.set_ylabel("Effective reproduction number - Rt")


def model_epidemic(data_file, mean_si, sd_si, window = 1, plot_start_day = 7,  uncertain_w = False,  plot_w = False, plot_incidence = False, plot_r_t = False, plot_surface = None, label = '', with_prediction = False):
    days, incidence_data = read_file(data_file)
    #incidence_data = [np.random.negative_binomial(i, 0.2) + i - 1 for i in incidence_data]
    w = infection_profile(mean_si, sd_si)
    T = days[-1]
    serial_interval_prob = []
    for t in range(T):
        serial_interval_prob.append(w(t))

    if plot_incidence:
        if not with_prediction:
            plot_incidence_data(range(T), incidence_data, label=label)
            return
        else:
            #plot_surface.scatter(range(1, T+ 1), incidence_data, zorder=20, marker="x", label="Actual Readings")
            #plot_surface.plot(range(1, T+1), incidence_data, zorder=0, label="test data")
            train_data = incidence_data[:-6]
            test_data = incidence_data[-6:]
            predicted = []
            estimates_of_R_t = []
            current_day = 0

            for t in range(len(train_data)):
                estimates_of_R_t.append(estimate_R_t(current_day, window, incidence_data=train_data, w=w))
                current_day += 1

            plot_surface.scatter(range(1, len(train_data) + 1), train_data, zorder = 10, label = "Reported cases", c='b')
            plot_surface.scatter(range(len(train_data), len(train_data) + len(test_data) + 1), [train_data[-1]] + test_data, zorder = 5, c='b')

            for t in range(len(test_data)):
                estimates_of_R_t.append(estimate_R_t(current_day, window, incidence_data=train_data,
                                                     w=w))  # getting r_t for next prediction
                predicted.append(predict_I_t(current_day, estimates_of_R_t[-1][-1], incidence_data = train_data, w=w)) #getting our prediction using last value of R_t
                train_data.append(test_data[t]) #adding actual value to train data
                current_day += 1

            estimates_of_R_t.append(estimate_R_t(current_day, window, incidence_data=train_data,
                                                 w=w))
            predicted.append(predict_I_t(current_day, estimates_of_R_t[-1][-1], incidence_data=train_data,
                                         w=w))  # getting our prediction using last value of R_t
            train_data.append(predicted[-1][-1])




            starts, ends, means = [],[], []

            for t in range(len(predicted)):
                means.append(predicted[t][-1])
                starts.append(predicted[t][0][0])
                ends.append(predicted[t][0][1])
            plot_surface.fill_between(range(len(train_data) - len(predicted) + 1, len(train_data) + 1), starts, ends, alpha = 0.1, color = "black", zorder = 10, label = "confidence interval of prediction")
            plot_surface.plot(range(len(train_data) - len(predicted), len(train_data) + 1), [train_data[-len(predicted)-1]] + means, 'r', zorder = 15, color = "black", label = "mean of prediction")
            plot_surface.set_xlabel("Days T")
            plot_surface.set_ylabel("No of New cases")



            means, start, end = [], [], []
            for i in range(T):
                m = estimates_of_R_t[i][-1]
                means.append(m)
                start.append(estimates_of_R_t[i][0][0])
                end.append(estimates_of_R_t[i][0][1])


    if not uncertain_w:
        if plot_w:
            plot_serial_interval(range(1, T + 1), serial_interval_prob, label = "mean = " + str(mean_si) + " SD = " + str(sd_si))

        elif plot_r_t:
            estimates_of_R_t = []
            for t in range(T):
                estimates_of_R_t.append(estimate_R_t(t, window, incidence_data=incidence_data, w=w))
            means, start, end = [], [], []
            for i in range(T):
                m = estimates_of_R_t[i][-1]
                means.append(m)
                start.append(estimates_of_R_t[i][0][0])
                end.append(estimates_of_R_t[i][0][1])

            plot_estimates_r_t(range(plot_start_day, T), starts=start[plot_start_day:], ends=end[plot_start_day:], means=means[plot_start_day:], label=label)
        return

    if uncertain_w:
        serial_interval_prob = []
        N = 1000
        means = truncated_normal(mean_si, lower= 3.7, upper= 6.0, sd = 1, n = N)
        sd = [max(means) + 1 for i in range(N)]
        for i in range(N):
            while sd[i] >= means[i]:

                sd[i] = truncated_normal(sd_si, lower=1.9, upper= 4.9, sd = 1)
        mean_of_means = round(np.mean(means), 2)
        mean_of_SD = round(np.mean(sd), 2)

        print(mean_of_means, mean_of_SD)
        w = infection_profile(mean_of_means, mean_of_SD)

        if plot_w:
            for t in range(T):
                serial_interval_prob.append(w(t))
            plot_serial_interval(range(1, T + 1), serial_interval_prob, label = "mean = " + str(mean_of_means) + " SD = " + str(mean_of_SD))

        elif plot_r_t:
            estimates_of_r_t = np.zeros((N, T))
            for n in range(N):
                print(n)
                w = infection_profile(means[n], sd[n])
                for t in range(T):
                    estimates_of_r_t[n, t] = estimate_R_t(t = t, pi = window, incidence_data = incidence_data, w = w, sample=True, n = 1)[0]

            starts_of_r_t, ends_of_r_t, means_of_r_t = [],[],[]
            for t in range(plot_start_day, T):
                data = estimates_of_r_t[:, t]
                mean, start, end = mean_confidence_interval(data)
                starts_of_r_t.append(start)
                ends_of_r_t.append(end)
                means_of_r_t.append(mean)

            plot_estimates_r_t(range(plot_start_day, T), starts=starts_of_r_t, ends=ends_of_r_t,
                               means= means_of_r_t, label=label)
        return
















    """
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


        plot_estimates_r_t(range(day, T), starts = start[day:], ends=end[day:], means= means[day:], label = label)
    """



names = ["Punjab", "Sindh", "GB", "ICT", "KPK", "AJK", "Balochistan", "Pakistan", "United Kingdom", "Italy", "Spain", "China"]
names = [i + ".txt" for i in names]

select = [i for i in names]



WINDOW = 2
MEAN_SERIAL_INTERVAL = 8.4
STD_SERIAL_INTERVAL = 3.8



for i in range(len(select)):
    fig, axs = plt.subplots(2, 1)
    name = select[i]
    plot_surface = axs[1]
    model_epidemic(name, plot_start_day= 20, window = WINDOW, mean_si= MEAN_SERIAL_INTERVAL, sd_si= STD_SERIAL_INTERVAL, plot_incidence =  True, plot_surface= plot_surface, with_prediction = True)
    start_date = datetime.date(2020, 2, 28)
    delta = datetime.timedelta(days=10)
    dates = []
    for i in range(51):
        dates.append(start_date)
        start_date += delta
    dates = [i.strftime("%d-%b") for i in dates]
    plot_surface.set_xticklabels(dates)
    #plot_surface.axhline(1, ls='--', label = '1')
    plot_surface.legend()
    plot_surface.grid()
    #plot_surface.set_ylim((0, 5))
    plot_surface = axs[0]
    model_epidemic(name, plot_start_day= 20, window = WINDOW, mean_si= MEAN_SERIAL_INTERVAL, sd_si= STD_SERIAL_INTERVAL, plot_r_t =  True, plot_surface= plot_surface)
    plot_surface.legend()
    plot_surface.set_title(name[:-4].capitalize() + " Data with window size = " + str(WINDOW))
    start_date = datetime.date(2020, 2, 28)
    delta = datetime.timedelta(days=5)
    dates = []
    for i in range(51):
        dates.append(start_date)
        start_date += delta
    dates = [i.strftime("%d-%b") for i in dates]
    plot_surface.set_xticklabels(dates[5:])
    plot_surface.grid()
    # plt.tight_layout()
    fig.set_size_inches(16, 9)
    plt.savefig("Predictions\\"+ name[:-4].capitalize() + "Predicted.pdf", bbox_inches  = 'tight', dpi = 100)
    #fig.autofmt_xdate()

    #plt.show()
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
