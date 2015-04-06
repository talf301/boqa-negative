import hpo

import logging
from omim import MIM


class Net:
    pass

class ItemNode:
    pass

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

