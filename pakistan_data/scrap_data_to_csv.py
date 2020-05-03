from bs4 import BeautifulSoup
import datetime
from requests import get
from pytz import timezone
import csv

my_date = datetime.datetime.now(timezone('Etc/GMT+5'))
current_date = my_date.date()
current_date_str = my_date.date().strftime('%d-%b-%y')






mapping = {"Sindh": 2, "Punjab": 1, "KP": 3, "Islamabad":0, "Balochistan": 4, "GB": 5, "AJK":6}

def scrap_data():
    """to scrap current data"""
    url = 'https://www.geo.tv/'
    response = get(url)

    html_soup = BeautifulSoup(response.text, 'html.parser')

    corona_banner = html_soup.find('div', class_='coronavirus_banner')
    corona_banner = corona_banner.find_all('ul')

    data = [0 for i in range(len(mapping))]

    for ul in (corona_banner):
        for li in ul.find_all('li'):
            name = (li.find("span", class_ = "vb_right_number").text)
            number = int(li.find("span", class_="vb_right_text").text)
            if name in mapping:
                data[mapping[name]] = number
    return data

def write_to_csv(data):
    file_name = 'pakistan_data.csv'
    with open(file_name, newline='', mode= 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        rows = []
        for row in spamreader:
            rows.append(row)

    last_reading = rows[-1]
    last_date = (datetime.datetime.strptime(last_reading[0], '%d-%b-%y'))
    last_date = last_date.date()

    if last_date == current_date:
        rows[-1] = [current_date.strftime('%d-%b-%y')] + data
    else:
        assert (abs(current_date - last_date).days) == 1
        rows.append( [current_date.strftime('%d-%b-%y')] + data)


    file_name = 'pakistan_data.csv'
    with open(file_name, newline='', mode='w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|')
        for row in rows:
            spamwriter.writerow(row)


if __name__ == "__main__":
    write_to_csv(scrap_data())


