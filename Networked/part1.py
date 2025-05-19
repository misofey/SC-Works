from model import A, B
import numpy as np
import control as ct

def Q1():
    poles = [-4.9, -5.1]
    K = ct.place(A, B, poles)


if __name__ == "__main__":
    Q1()
