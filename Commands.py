p20.transfer(1, R1['WELL1'], MISC['OUT'])
p20.transfer(1, R1['WELL2'], MISC['OUT'])
p20.transfer(1, R2['WELL3'], MISC['OUT'])
p20.transfer(1, R2['WELL4'], MISC['OUT'])
flag10mer = True
if flag10mer:
    p20.transfer(1, MISC['A1'], MISC['OUT']) # Backbone
    p20.transfer(1, MISC['A2'], MISC['OUT']) # Buffer
    p20.transfer(3.25, MISC['A3'], MISC['OUT']) # Water
else:
    p20.transfer(1, R3['WELL5'], MISC['OUT'])
    p20.transfer(1, R3['WELL6'], MISC['OUT'])
    p20.transfer(1, MISC['A1'], MISC['OUT']) # Backbone
    p20.transfer(1, MISC['A2'], MISC['OUT']) # Buffer
    p20.transfer(1.25, MISC['A3'], MISC['OUT']) # Water