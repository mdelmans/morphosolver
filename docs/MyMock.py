import sys

class MockBase:
	def __init__(self, *args, **kwargs):
		pass

class MockType(type):
	def __new__(meta, name, parent = None, *args, **kwargs):
		mockType = super().__new__(meta, name, (MockBase,), {})

		if parent != name:
			mockType.__module__ = parent
		else:
			mockType.__module__ = None

		return mockType

	def __iter__(self):
		return self

	def __next__(self):
		raise StopIteration

	def __getitem__(self, indices):
		raise IndexError

	def __init__ (self, name, parent = None, *args, **kwargs):
		self._name = name
		self._parent = parent

	# def __call__(self, *args, **kwargs):
	# 	mockType = type(self._name, (), {})
	# 	mockType.__module__ = self._parent
	# 	return mockType()

	def __getattr__(self, key):
		if key == '__all__':
			return []

		if not key.startswith('_'):
			mockModule = self._name
			print ("Parent: ", self._parent, " Name: ", self._name)
			if self._parent:
				mockModule = self._parent + '.' + mockModule

			mockType = MockType(key, mockModule)
			return mockType
		else:
			return ''

sys.modules['foo'] = MockType('foo') 

# print(i(11))


import foo
from foo import bar

print( dir(foo) )

# # from foo import *

# # from foo import bar

# # from foo import bar


# # # m.__module__ = 'parent'

import pickle

# class bar:
# 	pass

# class myclass(foo.bar.bax):
# 	pass

# print("Foobar: ", foo.bar())

pickle.dumps(bar)
