#
# This file is part of GA (collection of classical Graph Algorithms).
# 
# GA is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any
# later version.
#
# GA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with Octave; see the file COPYING.  If not
# see <http://www.gnu.org/licenses/>.
# 
# Copyright (C) 2009-2014 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#

import logging
import sys
from igraph import *

class FileWriter:
        def __init__(self, filename):
                self.f = open(filename, "w") 

        def write(self, stack):
                self.f.write(','.join(map(str, [i for i in stack])) + '\n')
        
        def close(self):
                self.f.close()

class MemoryWriter:
        def __init__(self):
                self.cycles = []

        def write(self,cycle):
                self.cycles.append(list(cycle))

class FilteredMemoryWriter:
        def __init__(self,l):
                self.cycles = []
                self.l = l

        def write(self,c):
                if len(c) == self.l:
                        self.cycles.append(list(c))

class Cycles:
  
        def __init__(self, g, writer, maxlength):
                """
                        g: directed graph to be analyzed
                        B: dictionary of parents for remembering uninteresting paths
                        blocked: array of blocked vertices
                        Ak: SCC with least vertex of subgraph induced by {s, ..., n} (relative vertex indices)          
                        sids: vertex indices in the original graph
                        stack: stack for storing the path being explored
                        s: index least vertex in the original graph
                        writer: file writer where to output the cycles found
                        maxlength: do not search beyond this limit
                """
                self.g = g
                self.B = {}
                self.blocked = {}
                self.Ak = None
                self.sids = []
                self.stack = []
                self.s = 0
                self.writer = writer
                self.maxlength = maxlength
                self.log = logging.getLogger('Cycles')
                
        def unblock(self, u):
                """ 
                        unblock a vertex u 
                """
                #self.log.info("unblock: " + str(u))
                self.blocked[u] = False
                for v in self.B[u]:
                        self.B[u].remove(v)
                        if self.blocked[v]:
                                self.unblock(v)

        def cycle(self, u):
                """ 
                        get all elementary cycles from vertex u (relative vertex index in subgraph K)
                """
                f = False
                nu = self.sids[u]
                self.stack.append(nu)
                self.blocked[nu] = True
                # NB: least vertex is 0 in Ak
                for v in self.Ak.successors(u):
                        nv = self.sids[v]
                        if v == 0:
                                self.writer.write(self.stack)   
                                f = True
                        elif len(self.stack) < self.maxlength: 
                                if not self.blocked[nv]:
                                        if self.cycle(v):
                                                f = True
                        else:
                                f = True
                if f:
                        self.unblock(nu)
                else:
                        for v in self.Ak.successors(u):
                                nv = self.sids[v]
                                if not (nu in self.B[nv]):
                                        self.B[nv].append(nu)
                self.stack.pop()
                return f

        def get_least_vertex_scc_vertices(self, sub):
                """ get the SCC with least vertex in the specified subgraph """
                membership = sub.components().membership
                for i in range(0, len(membership)):
                        if membership.count(membership[i]) > 1:
                                return [j+self.s for j in range(len(membership)) if membership[j] == membership[i]]                     
        
        def cycles(self):
                """ find all cycles in graph g, up to the specified max length """
                while True:
                        # subgraph induced by {s, ..., n}
                        sub = self.g.subgraph(range(self.s, self.g.vcount()))
                        # get vertices of SCC with least vertex 
                        cids = self.get_least_vertex_scc_vertices(sub)
                        if cids:
                                self.s = cids[0]
                                self.sids = cids
                                self.Ak = self.g.subgraph(cids)
                                self.log.info("s: " + str(self.s))
                                for i in cids:
                                        self.blocked[i] = False
                                        self.B[i] = []
                                self.cycle(0)
                                self.s = self.s + 1
                        else:
                                break
                                
if __name__ == "__main__":
        
        g = Graph.Full(7, "directed")
        
        writer = FileWriter("cycles.txt")
        c = Cycles(g, writer, 3)

        c.log.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        #ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        c.log.addHandler(ch)
        
        c.cycles()
        writer.close()
