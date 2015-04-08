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
        hpo = HPO(hpo_file)

        # Initialize empty stuff
        self.hid_dict = {}
        self.query_dict = {}
        self.hids = []
        self.items = []
        self.quers = []

        # Initialize hidden/query nodes
        for hp in hpo:
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
            item.fix_diseases(self.hid_dict, hpo)

if __name__ == '__main__':
    net = Net('./hp.obo', './phenotype_annotation.tab')
    print "Finished!"
    print len(net.hids)
    print len(net.items)
    print len(net.quers)




