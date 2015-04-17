import json

__author__ = 'Tal Friedman (talf301@gmail.com)'

if __name__ == '__main__':
    json_file = open('./phenotips_2015-02-08_18-45.json', 'r')
    data = json.load(json_file)
    diagnosed_patients = [d for d in data if 'disorders' in d.keys()]
    for i, d in enumerate(diagnosed_patients):
        try:
            disorder = d['disorders'][0]['id']
        except KeyError:
            continue
        disorder = disorder.split(':')[1]

        pos_phenos = [f['id'] for f in d.get('features', []) if f['observed'] == 'yes']
        # Just in case, not using this atm
        neg_phenos = [f['id'] for f in d.get('features', []) if f['observed'] == 'no']
        if pos_phenos:
            write_file = open('./phenotips_patients/' + disorder + '_' + str(i) + '_hpo.txt', 'w')
            write_file.write(','.join(pos_phenos) + '\n')
