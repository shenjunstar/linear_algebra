from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 6

class Vector(object):

	CANNOT_NORMALNIZE_ZERO_VECTOR_MSG = 'Cannot normalnize the zero vector'
	NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component'

	"""docstring for Vector"""
	def __init__(self, coordinates):
		super(Vector, self).__init__()
		try:
			if not coordinates:
				raise ValueError
			self.coordinates = tuple([Decimal(x) for x in coordinates])
			self.dimension = len(self.coordinates)
		except ValueError as e:
			raise ValueError('The coordinates must be nonempty')
		except TypeError as e:
			raise TypeError('The coordinates must be iterable')

	def  __str__(self):
		return 'Vector:{}'.format(self.coordinates)

	def  __eq__(self,v):
		return self.coordinates == v.coordinates

	#两个向量相加
	def  plus(self,v):
		new_coordinates = [x+y for x,y in zip(self.coordinates,v.coordinates)]
		return Vector(new_coordinates)
	
	#两个向量相减
	def minus(self,v):
		new_coordinates = [x-y for x,y in zip(self.coordinates,v.coordinates)]
		return Vector(new_coordinates)

	#常量与向量想乘
	def times_scalars(self,c):
		new_coordinates = [c*x for x in self.coordinates]
		return Vector(new_coordinates)

	#向量的长度
	def magnitude(self):
		new_coordinates = [x**Decimal('2') for x in self.coordinates]
		return sum(new_coordinates).sqrt()

	#求给定向量的单位向量
	def normalnized(self):
		try:
			magnitude = self.magnitude()
			return self.times_scalars(Decimal('1.0')/magnitude)
		except ZeroDivisionError as e:
			raise Exception('magnitude mustn\'t be zero')

	#向量的点集
	def dot(self,v):
		return sum([Decimal(x)*Decimal(y) for x,y in zip(self.coordinates,v.coordinates)])

	#两个向量的角度（弧度）
	def angle_with(self,v,in_degrees=False):
		try:
			u1 = self.normalnized()
			u2 = v.normalnized()
			angle_in_radians = acos(u1.dot(u2))
			if in_degrees:
				degrees_per_radian = 180./ pi
				return angle_in_radians*degrees_per_radian 
			else:
				return angle_in_radians
		except Exception as e:
			raise e

	#两个向量是否正交
	def is_orthogonal_to(self,v,tolerance=1e-10):
		return abs(self.dot(v)) < tolerance

	#两个向量是否平行
	def is_parallel_to(self,v):
		return (self.is_zero() or v.is_zero() or self.angle_with(v) ==0 or self.angle_with(v)==pi )

	#判断向量是否为0向量
	def is_zero(self,tolerance=1e-10):
		return self.magnitude() < tolerance

	#求向量投影（v在basis上的投影向量）
	def component_parallel_to(self,basis):
		try:
			u = basis.normalnized()
			weight = self.dot(u)
			return u.times_scalars(weight)
		except Exception as e:
			if str(e) == self.CANNOT_NORMALNIZE_ZERO_VECTOR_MSG:
				raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
			else:
				raise e

	#求向量投影的法向量（v在basis上的投影向量）
	def component_orthogonal_to(self,basis):
		try:
			projection = self.component_parallel_to(basis)
			return self.minus(projection)
		except Exception as e:
			if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
				raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
			else:
				raise e

	def cross(self,v):
		try:
			x_1,y_1,z_1 = self.coordinates
			x_2,y_2,z_2 = v.coordinates
			new_coordinates = [ y_1*z_2 - y_2*z_1 ,
								-(x_1*z_2 - x_2*z_1) ,
								x_1*y_2 - x_2*y_1 ]
			return Vector(new_coordinates)
		except ValueError as e:
			msg = str(e)
			if msg == 'need more than 2 values to unpack':
				self_embedded_in_R3 = Vector(self.coordinates + ('0',))
				v_embedded_in_R3 = Vector(v.coordinates + ('0',))
				return self_embedded_in_R3.cross(v_embedded_in_R3)
			elif (msg == 'too many values to unpack' or 
				msg == 'need more than one value to unpack'):
				raise Exception('need 2 or 3 values to unpack')
			else :
				raise Exception(e)

	def area_of_triangle_with(self,v):
		return self.area_of_parallelogram_with(v)/Decimal('2.0')

	def area_of_parallelogram_with(self, v):
		cross_product = self.cross(v)
		return cross_product.magnitude()
