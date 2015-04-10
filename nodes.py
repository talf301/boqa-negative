from omim import MIM, Disease
from hpo import HPO, HP
from bitarray import bitarray

__author__ = 'Tal Friedman (talf301@gmail.com)'

class ItemNode:
    """
    A single item node in the net

    Attributes:
    disease
    hidden nodes to which it is annotated w/ freq
    """

    def __init__(self, disease):
        """
        Create a new item node from a disease
        :param disease: Required omim disease object
        """
        self.disease = disease
        self.hids = {}

    def fix_diseases(self, hid_dict, hpo):
        """
        Update the dictionary of hidden units correctly
        """
        for hp_term, freq in self.disease.phenotype_freqs.items():
            # hp_term is a string, need to get the hp object and then hidden node
            try:
                hid_node = hid_dict[hpo[hp_term]]
                self.hids[hid_node] = freq
            except KeyError:
                pass

    def get_marginal_no_freq(self, hids, alpha, beta, interested_quer):
        """
        Compute the marginal for this disease ignoring frequency annotations
        """
        # hidden nodes actually annotated
        annot = bitarray(len(hids))
        annot.setall(False)

        for hid_node in self.hids.keys():
            annot[hid_node.index] = True
            for anc in hid_node.ancestors:
                annot[anc.index] = True
        #print len(annot) - len(annotated)
        # Do actual computation
        return self._compute_marginal(annot, hids, alpha, beta, interested_quer)


    def _compute_marginal(self, annotated, hids, alpha, beta, interested_quer):
        """
        Given an actual set of annotated hidden units (implicit and explicit), do computation
        :param annotated: bitarray of whether or not each unit is annotated
        :param hids: list of all hidden units in the net
        :param alpha: Alpha from model, false positive
        :param beta: Beta from model, false negative
        :param interested_quer: The terms which have all parents active in the query
        """
        m001 = 0
        m011 = 0
        m101 = 0
        m111 = 0
        for quer in interested_quer:
            # Only interested in units where visible parents are on
            hid = quer.hid
            if annotated[hid.index]:
                if hid.query.state == 1:
                    m111 += 1
                else:
                    m011 += 1
            else:
                if hid.query.state == 1:
                    m101 += 1
                else:
                    m001 += 1

        return pow(beta, m011) * pow(1-beta, m111) * pow(1-alpha, m001) * pow(alpha, m101)




class OntologyNode:
    """
    A single hidden or query node in the net

    Attributes:
    HP object
    parents
    children
    ancestors
    state
    """

    def __init__(self, hp, index):
        """
        Create a new Query/hidden node from an HP object
        :param hp: The hpo node that this is representing
        :param index: the index whih this node is given, used for "hashing" purposes
        """
        self.hp = hp
        self.parents = set()
        self.children = set()
        self.ancestors = set()
        self.state = None
        self.index = index

    def fix_parents_children(self, node_dict):
        """
        Assign the parent and children sets to be the set of OntologyNodes corresponding
        to the parent and children sets of the HP node

        :param node_dict: {hp -> OnotologyNode} for the appropriate layer
        """
        # Parents
        for term in self.hp.parents:
            self.parents.add(node_dict[term])

        # Children
        for term in self.hp.children:
            self.children.add(node_dict[term])

        # Ancestors
        for term in self.hp.ancestors():
            self.ancestors.add(node_dict[term])


class QueryNode(OntologyNode):
    """
    A single query node in the net

    Additional attributes:
    corresponding hidden node
    """

    def __init__(self, hp, index):
        """
        Create a new query node from an HP object
        :param hp: The hpo node being represented
        """
        self.hid = None
        self.parents_on = None
        OntologyNode.__init__(self, hp, index)

    def update_parent_status(self):
        self.parents_on = all(s.state == 1 for s in self.parents)


class HiddenNode(OntologyNode):
    """
    A single hidden node in the net
    Additional attributes:
    corresponding query node
    """

    def __init__(self, hp, index):
        """
        Create a new hidden node from an HP object
        :param hp:  The hpo node being represented
        :return:
        """
        self.query = None
        OntologyNode.__init__(self, hp, index)

