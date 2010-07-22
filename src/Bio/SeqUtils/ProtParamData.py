# This module contains indices to be used with ProtParam

# Kyte & Doolittle index of hydrophobicity
kd = {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
      'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
      'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
      'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

# Flexibility
# Normalized flexibility parameters (B-values), average (Vihinen et al., 1994)
Flex= {'A': 0.984, 'C': 0.906, 'E': 1.094, 'D': 1.068,
       'G': 1.031, 'F': 0.915, 'I': 0.927, 'H': 0.950,
       'K': 1.102, 'M': 0.952, 'L': 0.935, 'N': 1.048,
       'Q': 1.037, 'P': 1.049, 'S': 1.046, 'R': 1.008,
       'T': 0.997, 'W': 0.904, 'V': 0.931, 'Y': 0.929}

# Hydrophilicity
# 1 Hopp & Wood
# Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981). 
hw = {'A':-0.5,'R':3.0, 'N':0.2, 'D':3.0, 'C':-1.0,
      'Q':0.2, 'E':3.0, 'G':0.0, 'H':-0.5,'I':-1.8,
      'L':-1.8,'K':3.0, 'M':-1.3,'F':-2.5,'P':0.0,
      'S':0.3, 'T':-0.4,'W':-3.4,'Y':-2.3,'V':-1.5 }

# Surface accessibility
# 1 Emini Surface fractional probability
em = {'A':0.815,'R':1.475,'N':1.296,'D':1.283,'C':0.394,
      'Q':1.348,'E':1.445,'G':0.714,'H':1.180,'I':0.603,
      'L':0.603,'K':1.545,'M':0.714,'F':0.695,'P':1.236,
      'S':1.115,'T':1.184,'W':0.808,'Y':1.089,'V':0.606 }

# 2 Janin Interior to surface transfer energy scale
ja = {'A': 0.28,'R':-1.14,'N':-0.55,'D':-0.52,'C': 0.97,
      'Q':-0.69,'E':-1.01,'G': 0.43,'H':-0.31,'I': 0.60,
      'L': 0.60,'K':-1.62,'M': 0.43,'F': 0.46,'P':-0.42,
      'S':-0.19,'T':-0.32,'W': 0.29,'Y':-0.15,'V': 0.60 }


# A two dimentional dictionary for calculating the instability index.
# Guruprasad K., Reddy B.V.B., Pandit M.W.    Protein Engineering 4:155-161(1990).
# It is based on dipeptide values therefore the vale for the dipeptide DG is DIWV['D']['G'].
DIWV = {'A': {'A': 1.0, 'C': 44.94, 'E': 1.0, 'D': -7.49,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': -7.49,
              'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0,
              'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 1.0},
        'C': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 20.26,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 33.60,
              'K': 1.0, 'M': 33.60, 'L': 20.26, 'N': 1.0,
              'Q': -6.54, 'P': 20.26, 'S': 1.0, 'R': 1.0,
              'T': 33.60, 'W': 24.68, 'V': -6.54, 'Y': 1.0},
        'E': {'A': 1.0, 'C': 44.94, 'E': 33.60, 'D': 20.26,
              'G': 1.0, 'F': 1.0, 'I': 20.26, 'H': -6.54,
              'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': 1.0,
              'T': 1.0, 'W': -14.03, 'V': 1.0, 'Y': 1.0},
        'D': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0,
              'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 1.0, 'P': 1.0, 'S': 20.26, 'R': -6.54,
              'T': -14.03, 'W': 1.0, 'V': 1.0, 'Y': 1.0},
        'G': {'A': -7.49, 'C': 1.0, 'E': -6.54, 'D': 1.0,
              'G': 13.34, 'F': 1.0, 'I': -7.49, 'H': 1.0,
              'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': -7.49,
              'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0,
              'T': -7.49, 'W': 13.34, 'V': 1.0, 'Y': -7.49},
        'F': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 13.34,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0,
              'K': -14.03, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0,
              'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 33.601},
        'I': {'A': 1.0, 'C': 1.0, 'E': 44.94, 'D': 1.0,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 13.34,
              'K': -7.49, 'M': 1.0, 'L': 20.26, 'N': 1.0,
              'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0,
              'T': 1.0, 'W': 1.0, 'V': -7.49, 'Y': 1.0},
        'H': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': -9.37, 'F': -9.37, 'I': 44.94, 'H': 1.0,
              'K': 24.68, 'M': 1.0, 'L': 1.0, 'N': 24.68,
              'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0,
              'T': -6.54, 'W': -1.88, 'V': 1.0, 'Y': 44.94},
        'K': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': -7.49, 'F': 1.0, 'I': -7.49, 'H': 1.0,
              'K': 1.0, 'M': 33.60, 'L': -7.49, 'N': 1.0,
              'Q': 24.64, 'P': -6.54, 'S': 1.0, 'R': 33.60,
              'T': 1.0, 'W': 1.0, 'V': -7.49, 'Y': 1.0},
        'M': {'A': 13.34, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 58.28,
              'K': 1.0, 'M': -1.88, 'L': 1.0, 'N': 1.0,
              'Q': -6.54, 'P': 44.94, 'S': 44.94, 'R': -6.54,
              'T': -1.88, 'W': 1.0, 'V': 1.0, 'Y': 24.68},
        'L': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0,
              'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 33.60, 'P': 20.26, 'S': 1.0, 'R': 20.26,
              'T': 1.0, 'W': 24.68, 'V': 1.0, 'Y': 1.0},
        'N': {'A': 1.0, 'C': -1.88, 'E': 1.0, 'D': 1.0,
              'G': -14.03, 'F': -14.03, 'I': 44.94, 'H': 1.0,
              'K': 24.68, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': -6.54, 'P': -1.88, 'S': 1.0, 'R': 1.0,
              'T': -7.49, 'W': -9.37, 'V': 1.0, 'Y': 1.0},
        'Q': {'A': 1.0, 'C': -6.54, 'E': 20.26, 'D': 20.26,
              'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0,
              'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 1.0,
              'T': 1.0, 'W': 1.0, 'V': -6.54, 'Y': -6.54},
        'P': {'A': 20.26, 'C': -6.54, 'E': 18.38, 'D': -6.54,
              'G': 1.0, 'F': 20.26, 'I': 1.0, 'H': 1.0,
              'K': 1.0, 'M': -6.54, 'L': 1.0, 'N': 1.0,
              'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': -6.54,
              'T': 1.0, 'W': -1.88, 'V': 20.26, 'Y': 1.0},
        'S': {'A': 1.0, 'C': 33.60, 'E': 20.26, 'D': 1.0,
              'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0,
              'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 20.26, 'P': 44.94, 'S': 20.26, 'R': 20.26,
              'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 1.0},
        'R': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 20.26,
              'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 13.34,
              'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 58.28,
              'T': 1.0, 'W': 58.28, 'V': 1.0, 'Y': -6.54},
        'T': {'A': 1.0, 'C': 1.0, 'E': 20.26, 'D': 1.0,
              'G': -7.49, 'F': 13.34, 'I': 1.0, 'H': 1.0,
              'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': -14.03,
              'Q': -6.54, 'P': 1.0, 'S': 1.0, 'R': 1.0,
              'T': 1.0, 'W': -14.03, 'V': 1.0, 'Y': 1.0},
        'W': {'A': -14.03, 'C': 1.0, 'E': 1.0, 'D': 1.0,
              'G': -9.37, 'F': 1.0, 'I': 1.0, 'H': 24.68,
              'K': 1.0, 'M': 24.68, 'L': 13.34, 'N': 13.34,
              'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0,
              'T': -14.03, 'W': 1.0, 'V': -7.49, 'Y': 1.0},
        'V': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': -14.03,
              'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 1.0,
              'K': -1.88, 'M': 1.0, 'L': 1.0, 'N': 1.0,
              'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0,
              'T': -7.49, 'W': 1.0, 'V': 1.0, 'Y': -6.54},
        'Y': {'A': 24.68, 'C': 1.0, 'E': -6.54, 'D': 24.68,
              'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 13.34,
              'K': 1.0, 'M': 44.94, 'L': 1.0, 'N': 1.0,
              'Q': 1.0, 'P': 13.34, 'S': 1.0, 'R': -15.91,
              'T': -7.49, 'W': -9.37, 'V': 1.0, 'Y': 13.34},
        }
