#Codigo baseado no grupo 1

from Uniprot_Parser import *
from Bio.SeqIO import UniprotIO
import urllib


ids=['Q5F9Z1',
'Q5F9Z0',
'Q5F9Y9',
'Q5F9Y6',
'Q5F9Y5',
'Q5F9Y4',
'Q5F9Y3',
'Q5F9Y2',
'Q5F9Y1',
'Q5F9Y0',
'Q5F9X9',
'Q5F9X8',
'Q5F9X7',
'Q5F9X6',
'Q5F9X5',
'Q5F9X4',
'Q5F9X3',
'Q5F9X2',
'Q5F9X1',
'Q5F9X0',
'Q5F9W9',
'Q5F9W8',
'Q5F9W7',
'Q5F9W6',
'Q5F9W5',
'Q5F9W4',
'Q5F9W3',
'Q5F9W2',
'Q5F9W1',
'Q5F9W0',
'Q5F9V9',
'Q5F9V7',
'Q5F9V6',
'Q5F9V5',
'Q5F9V4',
'Q5F9V3',
'Q5F9V2',
'Q5F9V1',
'Q5F9V0',
'Q5F9U8',
'Q5F9U7',
'Q5F9U6',
'Q5F9U5',
'Q5F9U4',
'Q5F9U3',
'Q5F9U2',
'Q5F9U1',
'Q5F9U0',
'Q5F9T9',
'Q5F9T8',
'Q5F9T7',
'Q5F9T6',
'Q5F9T5',
'Q5F9T4',
'Q5F9T3',
'Q5F9T2',
'Q5F9T1',
'Q5F9T0',
'Q5F9S9',
'Q5F9S8',
'Q5F9S7',
'Q5F9S6',
'Q5F9S5',
'Q5F9S4',
'Q5F9S3',
'Q5F9S2',
'Q5F9S1',
'Q5F9S0',
'Q5F9R9',
'Q5F9R8',
'Q5F9R7',
'Q5F9R6',
'Q5F9R5',
'Q5F9R4',
'Q5F9R3',
'Q5F9R2',
'Q5F9R0',
'Q5F9Q9',
'Q5F9Q8',
'Q5F9Q7',
'Q5F9Q6',
'Q5F9Q5',
'Q5F9Q4',
'Q5F9Q3',
'Q5F9Q2',
'Q5F9Q1',
'Q5F9Q0',
'Q5F9P9',
'Q5F9P8',
'Q5F9P7',
'Q5F9P6',
'Q5F9P5',
'Q5F9P4',
'Q5F9P3',
'Q5F9P2',
'Q5F9P1',
'Q5F9P0',
'Q5F9N9',
'Q5F9N8',
'Q5F9N7',
'Q5F9N6',
'Q5F9N5',
'Q5F9N4',
'Q5F9N3',
'Q5F9N2',
'Q5F9N1',
'Q5F9N0',
'Q5F9M9',
'Q5F9M8',
'Q5F9M7',
'Q5F9M6',
'Q5F9M3',
'Q5F9M2',
'Q5F9M1',
'Q5F9M0',
'Q5F9L9',
'Q5F9L8',
'Q5F9L7',
'Q5F9L6',
'Q5F9L5',
'Q5F9L4',
'Q5F9L3',
'Q5F9L2',
'Q5F9L1',
'Q5F9L0',
'Q5F9K9',
'Q5F9K8',
'Q5F9K7',
'Q5F9K6',
'Q5F9K5',
'Q5F9K4',
'Q5F9K3',
'Q5F9K2',
'Q5F9K1',
'Q5F9K0',
'Q5F9J9',
'Q5F9J8',
'Q5F9J6',
'Q5F9J5',
'Q5F9J4',
'Q5F9J3',
'Q5F9J2',
'Q5F9J1',
'Q5F9J0',
'Q5F9I9',
'Q5F9I8',
'Q5F9I7',
'Q5F9I6',
'Q5F9I5',
'Q5F9I4',
'Q5F9I3',
'Q5F9I2',
'Q5F9I1',
'Q5F9I0',
'Q5F9H9',
'Q5F9H8',
'Q5F9H7',
'Q5F9H6',
'Q5F9H5',
'Q5F9H4',
'Q5F9H3',
'O87408',
'O87407',
'O87406',
'Q5F9G8',
'Q5F9G7',
'Q5F9G6',
'Q5F9G5',
'Q5F9G4',
'Q5F9G3',
'Q5F9G2',
'Q5F9G1',
'Q5F9G0',
'Q5F9F9',
'Q5F9F8',
'Q5F9F7',
'Q5F9F6',
'Q5F9F5',
'Q5F9F4',
'Q5F9F3',
'Q5F9F2',
'Q5F9F1',
'Q5F9F0',
'Q5F9E9',
'Q5F9E8',
'Q5F9E7',
'Q5F9E6',
'Q5F9E5',
'Q5F9E4',
'Q5F9E3',
'Q5F9E2',
'Q5F9E1',
'Q5F9E0',
'Q5F9D9',
'Q5F9D8',
'Q5F9D7',
'Q5F9D6',
'Q5F9D5',
'Q5F9D4',
'Q5F7P1',
'Q5F9D2',
'Q5F7P3',
'Q5F7P4',
'Q5F7P5',
'Q5F9C7',
'Q5F7P8',
'Q5F7P9',
'Q5F9C4',
'Q5F9C3',
'Q5F9C1',
'Q5F9C0',
'Q5F7Q6',
'Q5F7Q7',
'Q5F7Q8',
'Q5F9B5',
'Q5F9B4']


    
def more_info_uniprot():
    handle = open("proteoma4.xml")
    records=UniprotIO.UniprotIterator(handle,return_raw_comments=True)
    refs=[]
    for record in records:
        for i in range(len(ids)):
            if record.id==ids[i]:
                refs.append([ids[i]]+record.dbxrefs)#GO´s
    handle.close()
    return refs        
    
#print(more_info_uniprot())    


def tab():
    refs=more_info_uniprot()
    lista=[]
    l=[]
    for i in range(len(refs)):
        l.append(len(refs[i]))
        lista.append([])
    m=max(l)
    for j in range(len(refs)):
        if len(refs[j])<m:
            for k in range(len(refs[j])):
                lista[j].append(refs[j][k])
            for h in range(len(refs[j]),m):
                lista[j].append("nao tem")
        elif len(refs[j])==m:
            for g in range(len(refs[j])):
                lista[j].append(refs[j][g])        
    return lista
    
    
#sorting by protein ID
def sorting(lista):
    mat=[]
    for i in range(len(ids)):
        for j in range(len(lista)):
            if ids[i]==lista[j][0]:
                mat.append(lista[j])
    return mat
    

def GO():
    l=[]
    lista=sorting(tab())
    for i in range(len(lista)):
        for j in range(len(lista[i])):
            if "GO:GO:" in lista[i][j]:
                if lista[i][0] not in l:
                    l.append(lista[i][0])
                    l.append(lista[i][j])
                else:
                    l.append(lista[i][j])
                
                
                
    return l

print(GO())


