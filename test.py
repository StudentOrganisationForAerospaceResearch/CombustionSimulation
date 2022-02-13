import json

d = {
 "Laptop": {
            "sony": 1,
            "apple": 2,
            "asus": 5,
          },
 "Camera": {
            "sony": 2,
            "sumsung": 1,
            "nikon" : 4,
           },
}
with open("my.json","w") as f:
    json.dump(d,f)