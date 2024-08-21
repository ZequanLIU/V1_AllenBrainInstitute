import json

with open('V1_data.json', mode='r') as f:
    data = json.load(f)

print(data.keys())
print(data['simData'].keys())
print(data['net']['cells'][1])
print(data['net']['params']['cellParams']['VisL1Htr3a_0']['secs'].keys())
