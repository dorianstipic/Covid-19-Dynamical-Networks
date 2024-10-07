import numpy as np
from random import random

def R0simul(k_trip,vel_grozda,v):
    it=10000
    r=0.1
    days=16
    duljina=np.zeros(it)
    br_susjeda=vel_grozda - 1

    broj_izlazaka_p=np.zeros(it)
    p1_trip=1
    p2_trip=0.5
    #k_trip=5

    # simuliramo odkad je 1 (i samo 1) infectious u grozdu. 
    # simulacija traje days dana sto je broj dana inkubacije. Nakon toga dolazi karantena i ne gledamo.
    for k in range(it):
        grozd=np.zeros(vel_grozda) # svi u grozdu susceptible
        grozd[0]=1 # ubacimo jednog bolesnika na pocetku
        
        dana_bolestan=np.zeros(vel_grozda)
        dana_bolestan[0]=5
        for j in range(0,days):
            pomocni=grozd # uvedemo pomocni da se ne dogodi a-b, b-c u istom danu.
            for i in range(0,len(pomocni)): # idemo po svim susjedima prvog oboljelog
                if pomocni[i]==0: # ako je vrh susceptible ...
                    for i1 in range(int(np.sum(dana_bolestan>0))): # svaki bolesnik u grozdu
                        pomocni[i]+=int(random()<=r) # susceptible vrh ima sansu r da se zarazi (racunato za svakog bolesnog susjeda posebno)
                    pomocni[i]=min(1,pomocni[i]) # stavljamo novooboljele na vrijednost 1 a susceptible ostavimo na 0.
                    dana_bolestan[i]=5*min(1,pomocni[i])
                    if pomocni[i]==1: # svim vrhovima koji su se zarazili u ovoj iteraciji racunamo koliko dana izlaze.
                        if v=="version1":
                            broj_izlazaka_p[k]=p1_trip*(broj_izlazaka_p[k]+days-j)
                        if v=="version2":
                            if j<=5: # nakon toga slijedi karantena za grozd p2.
                                broj_izlazaka_p[k]=p1_trip*(broj_izlazaka_p[k]+5-j)+p1_trip*p2_trip*(broj_izlazaka_p[k]+days-5)
                                # unutar prvih 5 dana djeluje samo p1, kasnije p1*p2.
                            else:
                                broj_izlazaka_p[k]=p1_trip*p2_trip*(broj_izlazaka_p[k]+days-j)
                else:
                    dana_bolestan[i]-=1
                    continue
            grozd=pomocni
            
        if v=="version1":
            broj_izlazaka_p[k]=p1_trip*(broj_izlazaka_p[k]+days) # prvi bolesnik u grozdu ce putovati "days" dana.
        if v=="version2":
            broj_izlazaka_p[k]=p1_trip*(broj_izlazaka_p[k]+5)+p1_trip*p2_trip*(broj_izlazaka_p[k]+days-5) # prvi bolesnik u grozdu ce putovati 5 dana.
            
        duljina[k]=np.sum(grozd)-1 # novozarazeni u grozdu bez onog prvog

    broj_zarazenih_susjeda=np.mean(duljina) 
    broj_zarazenih_na_putovanju=np.mean(broj_izlazaka_p)*k_trip*r
    ukupan_broj_zarazenih=broj_zarazenih_susjeda+broj_zarazenih_na_putovanju
    R0=ukupan_broj_zarazenih/vel_grozda

    print("velicina grozda = ",vel_grozda, "\nbroj susjeda = ",br_susjeda,
          "\nprosjecan broj zarazenih susjeda = ",broj_zarazenih_susjeda,
          "\nbroj zarazenih na putovanjima = ", broj_zarazenih_na_putovanju,
         "\nUkupan broj zarazenih od grozdovih sudionika dok se grozd ne karantira = ",ukupan_broj_zarazenih,
         "\nR0 = ", R0)
    return R0
