import numpy as np
import cvxpy as cp
from itertools import product

KB = 1000#Bytes
MB = 1000*KB

PKT = 1500*KB

class LinkBase:
    LinkRate = 0
    AC0 = []
    AC1 = []
    AC2 = []
    AC3 = []

    @classmethod
    def iter(cls):
        return (cls.AC0, cls.AC1, cls.AC2, cls.AC3)

    @classmethod
    @property
    def name(cls):
        return cls.__name__
    pass

class AppBase:
    def __init__(self):
        self.variable = cp.Variable()
    pass

class RTApp(AppBase):
    def __init__(self, arrival, max_qos, weight=1.0):
        super().__init__()
        self.arrival = arrival
        self.max_qos = max_qos
        self.weight = weight
    @property
    def utility(self):
        return self.arrival*PKT
    @property
    def constraints(self) -> list:
        return [
            self.variable >= self.arrival*PKT,
            self.calc_qos() <= self.max_qos
        ]
    
    def calc_qos(self):
        return 0
    pass
# cp.multiply, cp.exp, cp.inv_pos, cp.sum
class DLApp(AppBase):
    def __init__(self, arrival, max_qos, weight=1.0):
        super().__init__()
        self.arrival = arrival
        self.max_qos = max_qos
        self.weight = weight
    @property
    def utility(self):
        return self.arrival*PKT
    @property
    def constraints(self) -> list:
        return [
            self.variable >= self.arrival*PKT,
            self.calc_qos() <= self.max_qos
        ]
    
    def calc_qos(self):
        return 0
    pass

class ThruApp(AppBase):
    def __init__(self, min_thru, size=np.Inf, weight=1.0):
        super().__init__()
        self.size = size
        self.min_thru = min_thru#Bps
        self.weight = weight
    @property
    def utility(self):
        return self.variable
    @property
    def constraints(self) -> list:
        return [
            self.variable >= self.min_thru
        ]
    def calc_qos(self):
        return - self.weight*self.variable
    pass

class ACBase:
    __slots__ = ['AIFSn', 'CWmin', 'TXOP']
    def __init__(self, phy='ht', mcs=54):
        self.mcs = mcs / 8 * MB #Bps
        match phy:
            case 'ht':  self.amsdu = 7935#Bytes
            case 'vht': self.amsdu = 11454#Bytes
            case _:     self.amsdu = 3839#Bytes
        pass
    
    @property
    def txop(self):
        if self.TXOP>0:
            return self.TXOP#ms
        else:
            return (self.amsdu/self.mcs)*1000#ms
    pass

class AC0(ACBase):
    AIFSn = 2
    CWmin = 3
    TXOP  = 1.504#ms

class AC1(ACBase):
    AIFSn = 2
    CWmin = 7
    TXOP  = 3.008#ms

class AC2(ACBase):
    AIFSn = 3
    CWmin = 15
    TXOP  = 0

class AC3(ACBase):
    AIFSn = 7
    CWmin = 15
    TXOP  = 0

def contend_portion(A:ACBase, B:ACBase) -> np.array:
    ##
    CWa, CWb     = A.CWmin, B.CWmin
    AIFSa, AIFSb = A.AIFSn, B.AIFSn
    ###
    vCWa = CWa + AIFSa
    vCWb = CWb + AIFSb
    P1 = np.zeros(( (CWa +AIFSa) * (CWb + AIFSb),(CWa +AIFSa) * (CWb + AIFSb) ))
    P2 = np.zeros(( (CWa +AIFSa) * (CWb + AIFSb),(CWa +AIFSa) * (CWb + AIFSb) ))

    ##
    for a,b in product( range(vCWa), range(vCWb) ):
        if a >= b :
            P1[a + b * vCWa , (a-b) ] = 1
        else:
            P1[a + b * vCWa, (b-a) * vCWa] = 1
    ###
    for a,b in product( range(CWa), range(CWb) ):
        if a == 0 and b != 0:
            for temp in range(CWa):
                P2[a + AIFSa + b * vCWa, temp + AIFSa + b * vCWa] = 1/CWa
        elif b == 0 and a != 0:
            for temp in range(CWb):
                P2[a + (b + AIFSb) * vCWa, a  + (temp + AIFSb) * vCWa] = 1/CWb
        elif a == 0 and b == 0:
            for temp_a in range( CWa ):
                for temp_b in range( CWb ):
                    P2[0 + 0 * CWa , (temp_a + AIFSa) + (temp_b + AIFSb) * CWa] = 1/CWa * 1/CWb
    
    ##
    result = np.dot(P1,P2)
    modified_matrix = result.transpose() - np.identity(vCWa* vCWb)
    modified_matrix = np.append(modified_matrix,np.ones((1, vCWa* vCWb)),axis = 0)
    ###
    vector_b = np.zeros(( vCWa* vCWb+1,1))
    vector_b[vCWa* vCWb,0] = 1
    C = np.linalg.lstsq(modified_matrix,vector_b, rcond=None)
    stationary_vector = C[0]
    ###
    probability = np.zeros(2)
    for a in range(vCWa):
        for b in range(vCWb):
            if a > b:
                probability[0] = probability[0] + stationary_vector[ a + b * vCWa,0]
            elif a < b:
                probability[1] = probability[1] + stationary_vector[a + b * vCWa,0]
    return probability

if __name__=='__main__':
    print( 'AC2:AC2 = ', contend_portion(AC2, AC2) )
    # print( 'AC1:AC2 = ', contend_portion(AC1, AC2) )
    # print( 'AC2:AC3 = ', contend_portion(AC2, AC3) )
    # print( 'AC1:AC3 = ', contend_portion(AC1, AC3) )