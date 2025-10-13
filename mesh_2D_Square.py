
from fem_student_Square import Node, Element

class Mesh2D:

    def __init__(self,nodes: "list",elements: "list"):
        """
        : param nodes : list of the initial mesh nodes
        : type : list [Node]
        : param nodes : list of the initial mesh elements
        : type : list [Element]
        """
        self.nodes=nodes
        self.elements=elements

        # param edges : dictionnay of the mesh edges
        # type : dictionnary {key : Integer, value : Edge}
        self.edges = {}

        # assign identifiers ids to mesh elements
        self.numerate_elements()

         # assign identifiers ids to mesh nodes
        self.numerate_nodes()

        # build mesh edges
        self.build_edges()



    def numerate_nodes(self):
        """
        Numerate nodes, the known nodes on Dirichlet boundary conditions are numerated last
        """
        node_id = 0
        # count free nodes first
        for node in self.nodes:
            if not (node.is_dirichlet):
                node.ID = node_id
                node_id += 1
        # count other nodes last
        for node in self.nodes:
            if (node.is_dirichlet):
                node.ID = node_id
                node_id += 1

    def numerate_elements(self):
        """
        Numerate elements. No specific order is applied
        """
        element_id = 0
        for element in self.elements:
            element.ID = element_id
            element_id += 1


    def refine(self, nb_dicho : "int"):
        for i in range(nb_dicho):
            self.split()



    def split(self):
        """
        refine the mesh by splitting 1 triangle into 4
        """
        # create new nodes
        new_nodes = {}
        for k in self.edges:
            edge = self.edges[k]
            new_node = Node((edge.node1.x + edge.node2.x) / 2., (edge.node1.y + edge.node2.y) / 2.)

            if edge.node1.is_dirichlet and edge.node2.is_dirichlet:
                if edge.node1.value == edge.node2.value:
                    new_node.set_dirichlet(edge.node1.value)
            new_nodes.update({k: new_node})
            self.nodes.append(new_node)

            if hasattr(Node, 'is_convection'):
                if edge.node1.is_convection() and edge.node2.is_convection() and edge.element2 is None:
                    new_node.set_convection((edge.node1.h + edge.node2.h) / 2, (edge.node1.t_inf + edge.node2.t_inf) / 2)
                elif edge.node1.is_convection() and edge.element2 is None:
                    new_node.set_convection(edge.node1.h, edge.node1.t_inf)
                elif edge.node2.is_convection() and edge.element2 is None:
                    new_node.set_convection(edge.node2.h, edge.node2.t_inf)

            if hasattr(Node, 'is_convection'):
                if edge.node1.is_source() and edge.node2.is_source():
                    new_node.set_source((edge.node1.q0 + edge.node2.q0) / 2)

        # make the new mesh
        new_mesh = []
        for element in self.elements:
            ext_nodes = element.nodes
            mid_nodes = []
            for i in range(len(element.nodes)):
                [n1, n2] = self.get_edge_nodes(element,i)
                mid_nodes.append(new_nodes.get(self.key(n1, n2)))
            # 4 new nodes
            e1 = Element([ext_nodes[0], mid_nodes[0], mid_nodes[2]],element.region)
            e2 = Element([mid_nodes[0], ext_nodes[1], mid_nodes[1]],element.region)
            e3 = Element([mid_nodes[1], ext_nodes[2], mid_nodes[2]],element.region)
            e4 = Element([mid_nodes[0], mid_nodes[1], mid_nodes[2]],element.region)
            new_mesh.append(e1)
            new_mesh.append(e2)
            new_mesh.append(e3)
            new_mesh.append(e4)

        # finalize
        self.elements = new_mesh
        self.numerate_nodes()
        self.numerate_elements()
        self.build_edges()

    def get_edge_nodes(self, element, i):
        """
        Return the nodes of an edge
        : param i: edge number
        : type : Integer
        : return: the 2 nodes of this edge
        """
        n1 = element.nodes[i]
        n2 = element.nodes[(i + 1) % len(element.nodes)]
        return [n1, n2]

    def build_edges(self):
        """
        Build or rebuild the edges. edges are unique
        """
        self.edges = {}
        for element in self.elements:
            for i in range(len(element.nodes)):
                [n1, n2] = self.get_edge_nodes(element, i)
                k = self.key(n1, n2)
                edge = self.edges.get(k)
                if edge is None:
                    edge = Edge(element, n1, n2)
                    self.edges.update({k: edge})
                else:
                    edge.set_second_element(element)

    @staticmethod
    def key(node1, node2):
        """
        return the key of an edge
        """
        if node1.ID < node2.ID:
            return str(node1.ID) + ' ' + str(node2.ID)
        else:
            return str(node2.ID) + ' ' + str(node1.ID)


class Edge:
    def __init__(self, element1, node1, node2):
        """
        Create a new edge
        : param element1 : first element sharing this edge
        : type : Element
        : param node1 : first node of the edge
        : type : Node
        : param node2 : second node of the edge
        : type : Node
        """
        self.element1 = element1
        if node2.ID < node1.ID:
            n=node2
            node2=node1
            node1=n

        self.node1 = node1
        self.node2 = node2
        self.element2 = None

    def set_second_element(self, element):
        """
        Set set second element shring the edge if any; Border edges do not have second element
        : param element : the second element
        : type : Element
        """
        self.element2 = element

    def __str__(self):
        """
        to string overloading
        :return: the string representation
        """
        return 'nodes[' + str(self.node1) + ' ' + str(self.node2) + '] elements [' + str(self.element1) + ' ' + str(
            self.element2) + ']'


