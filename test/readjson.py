import json



# read meshinfo
with open('setting.json') as json_file:
    meshinfo = json.load(json_file)
json_file.close()
print(meshinfo)