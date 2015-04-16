from hpo import HP, HPO
from omim import Disease, MIM
from nodes import HiddenNode, ItemNode, QueryNode
from bitarray import bitarray
import cProfile
import random

__author__ = 'Tal Friedman (talf301@gmail.com)'


class Net:
    def __init__(self, hpo_file, omim_file, neg_omim_file, ic_file=None):
        """
        Construct the layers of the net and then "attach" them
        :param hpo_file: path for hp.obo file
        :param omim_file: path for phenotype_annotation.tab file
        :param neg_omim_file: path for negative_phenotype_annotation.tab file
        :param ic_file: path for file with information contents
        :return:
        """

        # Create omim and hpo stuff
        mim = MIM(omim_file)
        diseases = [d for d in mim.diseases if d.db == 'OMIM']

        # Deal with loading negatives
        dis_dict = {dis.id: dis for dis in diseases}
        neg_omim = MIM(neg_omim_file)
        print len(neg_omim.diseases)
        for dis in neg_omim:
            try:
                dis_dict[dis.id].neg_pheno = dis.phenotype_freqs.keys()
            except KeyError:
                continue

        self.hpo = HPO(hpo_file)
        self.hpo.filter_to_descendants('HP:0000118')
        # Initialize empty stuff
        self.hid_dict = {}
        self.query_dict = {}
        self.hids = []
        self.items = []
        self.quers = []

        # Keep track of indices
        curr_ind = 0

        # Initialize hidden/query nodes
        for hp in self.hpo:
            # Create hidden/query node
            hid = HiddenNode(hp, curr_ind)
            quer = QueryNode(hp, curr_ind)
            curr_ind += 1
            # Add it to our {hp -> node} dicts
            self.hid_dict[hp] = hid
            self.query_dict[hp] = quer
            # Add it to our list of nodes
            self.hids.append(hid)
            self.quers.append(quer)

        # Dicts should line up
        assert self.hid_dict.keys() == self.query_dict.keys()

        # Initialize item nodes
        for dis in diseases:
            item = ItemNode(dis)
            self.items.append(item)

        # Fix hidden-query connections
        for hp in self.hid_dict.keys():
            hid = self.hid_dict[hp]
            quer = self.query_dict[hp]
            hid.query = quer
            quer.hid = hid
            hid.fix_parents_children(self.hid_dict)
            quer.fix_parents_children(self.query_dict)

        # Fix item-hidden connections
        for item in self.items:
            item.fix_diseases(self.hid_dict, self.hpo)

        # Load information content stuff if available
        self.ic_dict = {}
        if ic_file:
            self.ic_dict = self.load_ic(ic_file)

    def set_query(self, terms):
        """
        Set the states of the query layer based on a list of HPO ids
        :param terms: list of strings representing hpo ids
        :return:
        """
        # Set all to off first
        for quer in self.quers:
            quer.state = 0

        # Gather all to be activated
        to_activate = set()
        for term in terms:
            quer = self.query_dict[self.hpo[term]]
            to_activate.add(quer)
            for par in quer.parents:
                to_activate.add(par)

        # Activate
        for quer in to_activate:
            quer.state = 1

        # Update parent activation status
        for quer in self.quers:
            quer.update_parent_status()

        # Update the set of nodes we are actually interested in for computation
        self.interest_quer = [quer for quer in self.quers if quer.parents_on]
        self.interest_quer.sort()


    def diagnose(self, type='no_freq', k=5, n_samples=1000, p=0.001, alpha=0.001, beta=0.1):
        """
        Do the basic diagnosis, computing marginals and sorting and then returning the top 20 possibilities
        :return: a list of the top 5 diseases with probabilities
        """
        if type == 'no_freq':
            dis_scores = [(dis.disease.id
                           , dis.get_marginal_no_freq(self.hids, 0.001, 0.1, self.interest_quer)) for dis in self.items]
        elif type == 'k_freq':
            dis_scores = [(dis.disease.id,
                           dis.get_marginal_k_freq(self.hids, 0.001, 0.1, self.interest_quer, k=k)) for dis in self.items]
        elif type == 'sample':
            dis_scores = [(dis.disease.id,
                           dis.get_marginal_sampling(self.hids, 0.001, 0.1, self.interest_quer, samples=n_samples)) for dis in self.items]
        elif type == 'sample_p':
            # Generate samples first
            samples = []
            for i in xrange(n_samples):
                annot = bitarray(len(self.hids))
                annot.setall(False)
                for hid_node in self.hids:
                    if random.uniform(0,1) < p:
                        annot = annot | hid_node.bitarr
                samples.append(annot)
            dis_scores = [(dis.disease.id,
                           dis.get_marginal_sampling_p(self.hids, 0.001, 0.1, self.interest_quer, samples, n_samples=n_samples)) for dis in self.items]
        elif type == 'sample_ic':
            # Generate samples first
            samples = []
            for i in xrange(n_samples):
                annot = bitarray(len(self.hids))
                annot.setall(False)
                for hid_node in self.hids:
                    try:
                        if random.uniform(0,1) < p * self.ic_dict[hid_node]:
                            annot = annot | hid_node.bitarr
                    except KeyError:
                        continue
                samples.append(annot)
            dis_scores = [(dis.disease.id,
                           dis.get_marginal_sampling_p(self.hids, 0.001, 0.1, self.interest_quer, samples, n_samples=n_samples)) for dis in self.items]


        den = sum(d[1] for d in dis_scores)
        dis_scores = [(d[0], d[1]/den) for d in dis_scores]
        dis_scores.sort(key=lambda x: x[1], reverse=True)
        return dis_scores


    def load_ic(self, filename): # Return dict of HP code to info content (float)
        ic_dict = {}
        with open(filename) as f:
            for line in f:
                info = line.strip().split()
                try:
                    ic_dict[self.hid_dict[self.hpo['HP:' + info[0]]]] = float(info[1])
                except KeyError:
                    continue
        return ic_dict


if __name__ == '__main__':
    hp = './hp.obo'
    omim = './phenotype_annotation.tab'
    neg_omim = './negative_phenotype_annotation.tab'
    #cProfile.run('Net(hp, omim)')
    net = Net('./hp.obo', './phenotype_annotation.tab', './negative_phenotype_annotation.tab', './ic_parent.txt')
    print "Finished!"
    print len(net.hids)
    print len(net.items)
    print len(net.quers)
    print len(net.ic_dict.keys())
    net.set_query(open("./First_3450_356_hpo.txt", 'r').readline().split(','))
    #print(net.items[0].get_marginal_no_freq(net.hids, 0.001, 0.1))
    type = 'sample_ic'
    cProfile.run('net.diagnose(type=type)')
    print net.diagnose(type=type)[:20]




