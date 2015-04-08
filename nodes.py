from omim import MIM, Disease
from hpo import HPO, HP

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
            hid_node = hid_dict[hpo[hp_term]]
            self.hids[hid_node] = freq


class OntologyNode:
    """
    A single hidden or query node in the net

    Attributes:
    HP object
    parents
    children
    state
    """

    def __init__(self, hp):
        """
        Create a new Query/hidden node from an HP object
        :param hp: The hpo node that this is representing

        """
        self.hp = hp
        self.parents = set()
        self.children = set()
        self.state = False

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


class QueryNode(OntologyNode):
    """
    A single query node in the net

    Additional attributes:
    corresponding hidden node
    """

    def __init__(self, hp):
        """
        Create a new query node from an HP object
        :param hp: The hpo node being represented
        """
        self.hidden = None
        OntologyNode.__init__(self, hp)


class HiddenNode(OntologyNode):
    """
    A single hidden node in the net
    Additional attributes:
    corresponding query node
    """

    def __init__(self, hp):
        """
        Create a new hidden node from an HP object
        :param hp:  The hpo node being represented
        :return:
        """
        self.query = None
        OntologyNode.__init__(self, hp)

