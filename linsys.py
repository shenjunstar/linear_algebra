from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        self[row1],self[row2] = self[row2],self[row1]


    def multiply_coefficient_and_row(self, coefficient, row):
        n = self[row].normal_vector
        k = self[row].constant_term

        new_normal_vector = n.times_scalars(coefficient)
        new_constant_term = k * coefficient

        self[row] = Plane(normal_vector=new_normal_vector,constant_term=new_constant_term)


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term

        new_normal_vector = n1.times_scalars(coefficient).plus(n2)
        new_constant_term = (k1 * coefficient) + k2

        self[row_to_be_added_to] = Plane(normal_vector=new_normal_vector,constant_term=new_constant_term)


    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = system.dimension

        j=0
        for i in range(num_equations):
            while j < num_variables:
                c = MyDecimal(system[i].normal_vector.coordinates[j])
                if c.is_near_zero():
                    swap_succeeded = system.swap_with_row_below_for_nonzero_coefficient_if_able(i,j)
                    if not swap_succeeded:
                        j += 1
                        continue

                system.clear_coeffient_below(i,j)
                j += 1
                break

        return system

    def swap_with_row_below_for_nonzero_coefficient_if_able(self,row,col):
        num_equations = len(self)

        for k in range(row+1,num_equations):
            coefficient = MyDecimal(self[k].normal_vector.coordinates[col])
            if not coefficient.is_near_zero():
                self.swap_rows(row,k)
                return True

        return False

    def clear_coeffient_below(self,row,col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector.coordinates[col])

        for k in range(row+1,num_equations):
            n = self[k].normal_vector
            gamma = n.coordinates[col]
            alpha = -gamma/beta
            self.add_multiple_times_row_to_row(alpha,row,k)

    def compute_rref(self):
        tf = self.compute_triangular_form()

        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()

        for i in range(num_equations)[::-1]:
            j = pivot_indices[i]
            if j<0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(i,j)
            tf.clear_coeffient_above(i,j)

        return tf

    def scale_row_to_make_coefficient_equal_one(self,row,col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n.coordinates[col]
        self.multiply_coefficient_and_row(beta,row)

    def clear_coeffient_above(self,row,col):
        for k in range(row)[::-1]:
            n = self[k].normal_vector
            alpha = -(n.coordinates[col])
            self.add_multiple_times_row_to_row(alpha,row,k)

    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_extract_solution()
        except Exception as e:
            raise e

    def do_gaussian_elimination_and_extract_solution(self):
        rref = self.compute_rref()

        rref.raise_exception_if_contradictory_equation()
        rref.raise_exception_if_too_many_pivots()

        num_variables = rref.dimension
        solution_coordinates = [rref.planes[i].constant_term for i in range(num_variables)]

        return Vector(solution_coordinates)

    def raise_exception_if_contradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception('No Solution')
                else:
                    raise e

    def raise_exception_if_too_many_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots < num_variables:
            raise Exception('Inf solutions')


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

class Parametrization(object):
    """docstring for Parametrization"""
    def __init__(self, basepoint, direction_vectors):
        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.direction_vectors.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension
        except AssertionError:
            raise Exception('The basepoint and direction vectors should all live in the same dimension')
        

p1 = Plane(normal_vector=Vector(['5.262','2.739','-9.878']), constant_term='-3.441')
p2 = Plane(normal_vector=Vector(['5.111','6.358','7.638']), constant_term='-2.152')
p3 = Plane(normal_vector=Vector(['2.016','-9.924','-1.367']), constant_term='-9.278')
p4 = Plane(normal_vector=Vector(['2.167','-13.543','-18.883']), constant_term='-10.567')
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_solution()
print(r)