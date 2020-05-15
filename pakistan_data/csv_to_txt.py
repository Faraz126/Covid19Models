import csv
import datetime

def csv_to_txt(csv_file):
    with open(csv_file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        numbers = {"Pakistan": []}
        prov_names = []
        start_date, end_date = str(), str()

        num = 0
        for row in spamreader:
            if num == 0:
                prov_names = row[1:]
                for i in prov_names:
                    numbers[i] = []
            else:
                data = row[1:] #first two index have date
                if not start_date:
                    start_date = row[0]
                if not data:
                    print(f"File {csv_file} ended")
                    return
                for i in range(len(prov_names)):
                    numbers[prov_names[i]].append(- sum(numbers[prov_names[i]]) + int(data[i]))
                numbers["Pakistan"].append(- sum(numbers["Pakistan"]) + sum([int(i) for i in data]))

            end_date = row[0]
            num += 1



    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write(f"{start_date},{end_date}\n")
        file.writelines(lines)


def world_data_to_txt(country_name):
    countires = ['China', 'US', 'United Kingdom', 'Italy', 'France', 'Germany', 'Spain', 'Iran']
    country_index=  countires.index(country_name)
    with open("worldData.csv", newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        spamreader.__next__() #skipping first line
        numbers = {country_name: []}
        epidemic_start = False
        start_date = str()
        end_date = str()

        for row in spamreader:
            if int(row[country_index + 1]) != 0 and not epidemic_start:
                start_date = row[0]
                num = 1
                epidemic_start = True
            if epidemic_start:
                #print(num, row[country_index + 1])
                numbers[country_name].append(- sum(numbers[country_name]) + int(row[country_index + 1]))
                num += 1
            end_date = row[0]

    #print(numbers)
    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write(f"{start_date},{end_date}\n")
        file.writelines(lines)

def mobility_data_to_txt(country_names):
    """
    :param country_names: list of country names to extract mobility data for.
    """
    assert type(country_names) == list
    countries_with_data = {i.lower():[] for i in country_names}
    with open("global_mobility.csv", newline='', encoding='utf-8') as csvfile:

        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        header = next(spamreader) #first line
        print(header)
        for row in spamreader:
            country_name = row[1].lower()
            if country_name in countries_with_data and not row[2] and not row[3]:
                date = datetime.datetime.strptime(row[4], '%Y-%m-%d')
                countries_with_data[country_name].append((date, row[5:]))

    for country in countries_with_data:
        data_for_country = sorted(countries_with_data[country])
        txt_file = open(f'{country}_mobility.txt', "w")
        txt_file.write(header[4] + "\t" + "\t".join(header[5:]) + "\n")
        for row in data_for_country:
            current_date = row[0].strftime('%d-%b-%y')
            readings = "\t".join(row[1])
            txt_file.write(current_date + "\t" + readings + "\n")
        txt_file.close()




if __name__ == "__main__":
    #csv_to_txt("pakistan_data.csv")
    #world_data_to_txt("China")
    #world_data_to_txt("Spain")
    #world_data_to_txt("United Kingdom")
    #world_data_to_txt("Italy")
    mobility_data_to_txt(["United Kingdom", "Italy", "Spain"])


