from hpo import HP, HPO
from omim import Disease, MIM
from nodes import HiddenNode, ItemNode, QueryNode
import logging



class Net:
    def __init__(self, hpo_file, omim_file):
        """
        Construct the layers of the net and then "attach" them
        :param hpo_file: path for hp.obo file
        :param omim_file: path for phenotype_annotation.tab file
        :return:
        """
        # Create omim and hpo stuff
        mim = MIM(omim_file)
        diseases = [d for d in mim.diseases if d.db == 'OMIM']
        self.hpo = HPO(hpo_file)
        self.hpo.filter_to_descendants('HP:0000118')
        # Initialize empty stuff
        self.hid_dict = {}
        self.query_dict = {}
        self.hids = []
        self.items = []
        self.quers = []

        # Initialize hidden/query nodes
        for hp in self.hpo:
            # Create hidden/query node
            hid = HiddenNode(hp)
            quer = QueryNode(hp)
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
            quer.hidden = hid
            hid.fix_parents_children(self.hid_dict)
            quer.fix_parents_children(self.query_dict)

        # Fix item-hidden connections
        for item in self.items:
            item.fix_diseases(self.hid_dict, self.hpo)

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


    def diagnose(self):
        """

        :return: a list of the top 5 diseases with probabilities
        """
        dis_scores = [(dis.disease.id
                       , dis.get_marginal_no_freq(self.hids, 0.001, 0.1)) for dis in self.items]
        dis_scores.sort(key=lambda x: x[1], reverse=True)
        return dis_scores[:20]


if __name__ == '__main__':
    net = Net('./hp.obo', './phenotype_annotation.tab')
    print "Finished!"
    print len(net.hids)
    print len(net.items)
    print len(net.quers)
    net.set_query(open("./First_3450_356_hpo.txt", 'r').readline().split(','))
    print(net.items[0].get_marginal_no_freq(net.hids, 0.001, 0.1))
    print net.diagnose()




