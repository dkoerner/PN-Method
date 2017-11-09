
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

import json


response = urlopen("http://127.0.0.1:5000/1.0/2.222/3.0")
parsed_json = json.load(response)











