import math

def normal_pdf(x, mean, std_dev):
    
    a = 1 / (std_dev * math.sqrt(2 * math.pi))
    b = math.exp(-((x - mean)**2) / (2 * std_dev**2))
    
    value = a * b
    
    return round(value, 5)
