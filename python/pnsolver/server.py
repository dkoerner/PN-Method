import os
#import redis
#import urlparse
from werkzeug.wrappers import Request, Response
from werkzeug.routing import Map, Rule
from werkzeug.exceptions import HTTPException, NotFound
from werkzeug.wsgi import SharedDataMiddleware
from werkzeug.utils import redirect
from jinja2 import Environment, FileSystemLoader

import json

import pnsolver
import numpy as np
import util

filename = "C:/projects/epfl/epfl17/python/pnsolver/test.pns"
pns = pnsolver.load_solution(filename)



class WSGITest(object):

	def __init__(self, config):
		self.url_map = Map([
							Rule('/<float:x>/<float:y>/<float:z>', endpoint='test'),
		])

	def dispatch_request(self, request):
		#return Response('Hello World!')
		adapter = self.url_map.bind_to_environ(request.environ)
		try:
			endpoint, values = adapter.match()
			return getattr(self, 'on_' + endpoint)(request, **values)
		except (HTTPException) as e:
			return e

	def wsgi_app(self, environ, start_response):
		request = Request(environ)
		response = self.dispatch_request(request)
		return response(environ, start_response)

	def __call__(self, environ, start_response):
		return self.wsgi_app(environ, start_response)

	def on_test(self, request, x, y, z):
		#print("!!!!!!!!!!!!!!!!{} {} {}".format(x, y, z))
		pWS = np.array([x, y, z])
		points = []
		#values = []
		#points.append((x, y, z))

		theta_list = np.arange(0.0, 1.0, 0.05)*np.pi
		phi_list = np.arange(0.0, 1.0, 0.05)*2.0*np.pi

		max_value = 0.0
		for j in range(phi_list.shape[0]):
			for i in range(theta_list.shape[0]):
				theta = theta_list[i]
				phi = phi_list[j]
				d = util.sphericalDirection(theta, phi)
				value = pns.eval(pWS, d)
				max_value = np.max([value, max_value])
				p = d*value
				points.append((p[0], p[1], p[2]))



		d = {
		'numPoints': len(points),
		'points': points,
		'max':max_value
		}

		json_string = json.dumps(d)
		#print("sending:\n{}".format(json_string))
		return Response( json_string, mimetype="application/text" )
	


def create_app(redis_host='localhost', redis_port=6379, with_static=True):
    app = WSGITest({
        'redis_host':       redis_host,
        'redis_port':       redis_port
    })
    if with_static:
        app.wsgi_app = SharedDataMiddleware(app.wsgi_app, {
            '/static':  os.path.join(os.path.dirname(__file__), 'static')
        })
    return app


if __name__ == '__main__':
    from werkzeug.serving import run_simple
    app = create_app()
    run_simple('127.0.0.1', 5000, app, use_debugger=True, use_reloader=True)