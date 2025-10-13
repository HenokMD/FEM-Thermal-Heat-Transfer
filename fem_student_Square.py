import numpy as np
from scipy import sparse
from scipy.sparse import linalg


class Node:

    def __init__(self,x : "float" ,y : "float", value : float =0.0,is_dirichlet:"bool"= False):
        self.x=x
        self.y=y
        self.value=0.0
        self.ID=None
        self.is_dirichlet=False   # state of the Dirichlet type boundary condition at the node, True : the value of the solution in the node known, otherwise : False
        self.is_convection=False  # state of the convection type boundary condition at the node
        #self.set_dirichlet1(1)    # This method is called for every node created, regardless of its position.
        self.h = 0.0                # param h : the heat convection factor
        self.t_inf = 0.0          # param t_inf : temperature at infinity
        self.is_source = False
        self.q0 = 0.0             # param q0 : thermal flux source


    def set_dirichlet(self, value):
        """
        Set the node as a Dirichlet boundary node with a known value
        :param value: the prescribed value
        """
        self.is_dirichlet = True
        self.value = value

    def set_value(self,val : "float"):
        self.value=val

    def getCoordinates(self):
        return [self.x,self.y]


    def set_convection(self, h, ambient):
        """
        Set the node as a Neumann boundary node
        :param h: the convection factor
        :param ambient: the ambient temperature
        """
        if self.y==0:
            self.is_convection = True
            self.h = h
            self.t_inf = ambient

    def set_source(self, q):
        """
        Set a source
        :param q: source
        :return: None
        """
        self.q0 = q
        self.is_source = True

class Region:
    def __init__(self, name, sigma, source):
        """
        Define a new region
        :param name name of the region
        :param sigma: conductivity
        :param source: source term
        """
        self.name = name
        self.sigma = sigma
        self.source = source

    def __str__(self):
        """
        to string overloading
        :return: the string representation
        """
        return 'Region ' + self.name, ' sigma ' + str(self.sigma) + ' source ' + str(self.source)


class Element:

    def __init__(self, nodes : "list", region : "Region"):
        self.nodes=nodes
        self.region=region
        self.lamda=region.sigma
        self.D=(nodes[2].getCoordinates()[1]-nodes[1].getCoordinates()[1])*(nodes[1].getCoordinates()[0]-nodes[0].getCoordinates()[0])-(nodes[2].getCoordinates()[0]-nodes[1].getCoordinates()[0])*(nodes[1].getCoordinates()[1]-nodes[0].getCoordinates()[1])
        self.S = self.D/2

    def grad_phi1x(self):
        return (self.nodes[1].getCoordinates()[1]-self.nodes[2].getCoordinates()[1])/self.D

    def grad_phi2x(self):
        return (self.nodes[2].getCoordinates()[1]-self.nodes[0].getCoordinates()[1])/self.D

    def grad_phi3x(self):
        return (self.nodes[0].getCoordinates()[1]-self.nodes[1].getCoordinates()[1])/self.D

    def grad_phi1y(self):
        return (self.nodes[2].getCoordinates()[0]-self.nodes[1].getCoordinates()[0])/self.D

    def grad_phi2y(self):
        return (self.nodes[0].getCoordinates()[0]-self.nodes[2].getCoordinates()[0])/self.D

    def grad_phi3y(self):
        return (self.nodes[1].getCoordinates()[0]-self.nodes[0].getCoordinates()[0])/self.D

    def m00(self):
        return self.S*self.lamda*((self.grad_phi1x())**2+(self.grad_phi1y())**2)

    def m11(self):
        return self.S*self.lamda*((self.grad_phi2x())**2+(self.grad_phi2y())**2)

    def m22(self):
        return self.S*self.lamda*((self.grad_phi3x())**2+(self.grad_phi3y())**2)

    def m01(self):
        return self.S*self.lamda*(self.grad_phi1x()*self.grad_phi2x()+self.grad_phi1y()*self.grad_phi2y())

    def m02(self):
        return self.S*self.lamda*(self.grad_phi1x()*self.grad_phi3x()+self.grad_phi1y()*self.grad_phi3y())

    def m12(self):
        return self.S*self.lamda*(self.grad_phi3x()*self.grad_phi2x()+self.grad_phi3y()*self.grad_phi2y())