class Monomial(FixedTermPolynomial, valid_term_counts=(0, 1)):
    """Implements a single-variable monomial. A single-term polynomial."""

    def __init__(self, coefficient=1, degree=1):
        """Initialize the following monomial: coefficient * x^(degree)."""
        if not isinstance(degree, int):
            raise ValueError("Monomial's degree should be a natural number.")
        if degree < 0:
            raise ValueError("Polynomials cannot have negative-degree terms.")
        self._degree = degree
        self._coeff = coefficient

    def _trim(self):
        """Trims self._vector to length. Keeps constant terms."""

    @property
    def terms(self):
        """Get the terms of self as a list of tuples in coeff, deg form.

        Terms are returned from largest degree to smallest degree, excluding
        any terms with a zero coefficient.
        """
        if self._coeff == 0:
            return [(0, 0)]
        if self._degree == -inf:
            return [(0, 0)]
        return [(self._coeff, self._degree)]

    @terms.setter
    def terms(self, terms):
        """Set the terms of self as a list of tuples in coeff, deg form."""
        if not terms:
            self._coeff = 0
        elif len(terms) == 1:
            self._coeff, self._degree = terms[0]
        else:
            terms = sorted([term for term in terms if term[0] != 0], key=lambda x: x[1])
            if terms[0][1] == terms[-1][1]:
                self._coeff = sum(term[0] for term in terms)
                self._degree = terms[0][1]
            else:
                curr_coeff, curr_deg = terms[0]
                termx = []
                for coeff, deg in terms[1:]:
                    if curr_deg == deg:
                        curr_coeff += coeff
                    elif curr_coeff != 0:
                        if termx:
                            raise TermError("terms has more than one non-zero term.")
                        termx.append((curr_coeff, curr_deg))
                        curr_coeff = coeff
                        curr_deg = deg
                if termx:
                    if curr_coeff:
                        raise TermError("terms has more than one non-zero term.")
                    self._coeff, self._degree = termx[0]
                self.coeff = curr_coeff
                self.degree = curr_deg

    @property
    def _vector(self):
        """Get _vector."""
        if self.degree == -inf:
            return [0]
        return [0] * self._degree + [self._coeff]

    @_vector.setter
    def _vector(self, _vector):
        """Set _vector."""
        max_deg = len(_vector) - 1
        is_set = False
        for index, coeff in enumerate(reversed(_vector)):
            if coeff != 0:
                if is_set:
                    raise TermError("_vector has > 1 non-zero term.")
                self._coeff = coeff
                self._degree = max_deg - index
                is_set = True
        if not is_set:
            self._coeff = 0

    @classmethod
    def zero_instance(cls):
        """Return the Monomial which is 0."""
        return Monomial(0, 0)

    @property
    def coefficient(self):
        """Return the coefficient of the monomial."""
        return self._coeff

    @coefficient.setter
    def coefficient(self, coeff):
        """Set the coefficient of the monomial."""
        self._coeff = coeff

    @property
    def degree(self):
        """Return the degree of the monomial."""
        if self._coeff == 0:
            self._degree = -inf
        elif self._degree == -inf:
            self._degree = 0

        return self._degree

    @degree.setter
    def degree(self, degree):
        """Set the degree of the monomial."""
        self._degree = degree

    @extract_polynomial
    def __mul__(self, other):
        """Return self * other.

        The class which is more permissive will be returned.
        """
        if isinstance(other, Monomial) and self and other:
            return Monomial(self.coefficient * other.coefficient,
                self.degree + other.degree)
        return super().__mul__(other)

    @extract_polynomial
    def __rmul__(self, other):
        """Return other * self.

        The class which is more permissive will be returned.
        """
        return self * other

    def __lt__(self, other):
        """Return self < other.

        Compares the degrees of the monomials and then, if
        they are equal, compares their coefficients.
        """
        if self.degree == other.degree:
            return self.a < other.a
        return self.degree < other.degree

    def __gt__(self, other):
        """Return self > other.

        Compares the degrees of the monomials and then, if
        they are equal, compares their coefficients.
        """
        if self.degree == other.degree:
            return self.a > other.a
        return self.degree > other.degree

    def __pow__(self, power, modulo=None):
        """Return self ** power or pow(self, other, modulo)."""
        result = deepcopy(self)
        result **= power

        return result % modulo if modulo is not None else result

    def __ipow__(self, other):
        """Return self **= power.

        Assumes self is mutable.
        Does not mutate in the case that self == 0 and other != 1.
        """
        if not isinstance(other, int):
            raise ValueError(
                "Can't call Monomial() **= x with a non-integer type."
            )

        if other < 0:
            raise ValueError(
                "Monomial can only be raised to a non-negative power."
            )

        if not self:
            if other != 0:
                return self
            terms = [(1, 0)]
        else:
            terms = [(self.coefficient ** other, self.degree * other)]

        # No native option exists to modify Monomial degree.
        return self.try_set_self(terms)

    def __lshift__(self, other):
        """Return self << other.

        Returns a Monomial that is self * x^other.
        """
        if other < 0:
            return self >> -other

        if not self:
            return self.zero_instance()

        return Monomial(self.coefficient, self.degree + other)

    def __ilshift__(self, other):
        """Return self <<= other.

        Returns a Monomial that is self * x^other. Does not
        guarantee the same type is returned.
        """
        if other < 0:
            self >>= -other
            return self

        if not self:
            return self

        return self.try_set_self([(self.coefficient, self.degree + other)])

    def __irshift__(self, other):
        """Return self >>= other."""
        if other < 0:
            self <<= -other
            return self

        if not self:
            return self

        if other > self.degree:
            return self.try_set_self([(0, 0)])
        return self.try_set_self([(self.coefficient, self.degree - other)])

    def __repr__(self):
        """Return repr(self)."""
        deg = max(0, self.degree)
        return "Monomial({0!r}, {1!r})".format(self.coefficient, deg)

    def __getitem__(self, degree):
        """Get the coefficient of the term with the given degree."""
        if isinstance(degree, slice):
            return self._vector[degree]
        if degree == self.degree:
            return self._coeff
        if degree > self.degree or degree < 0:
            raise IndexError("Attempt to get coefficient of term with \
degree {0} of a {1}-degree monomial".format(degree, self.degree))
        return 0

    def __setitem__(self, degree, new_value):
        """Set the coefficient of the term with the given degree."""
        if isinstance(degree, slice):
            self._vector[degree] = new_value
        elif degree == self.degree:
            self.coefficient = new_value
        else:
            raise IndexError("Can not set more than 1 term on Monomial.")


class Constant(FixedDegreePolynomial, Monomial, valid_degrees=(0, -inf)):
    """Implements constants as monomials of degree 0."""

    def __init__(self, const=1):
        """Initialize the constant with value const."""
        Monomial.__init__(self, const, 0)

    @classmethod
    def zero_instance(cls):
        """Return the constant which is 0."""
        return Constant(0)

    def __eq__(self, other):
        """Return self == other."""
        if other == self.const:
            return True
        return super().__eq__(other)

    @property
    def const(self):
        """Return the constant term."""
        return self.coefficient

    @const.setter
    def const(self, val):
        """Set the constant term."""
        self._coeff = val

    @extract_polynomial
    def __mul__(self, other):
        """Return self * other."""
        if isinstance(other, Constant):
            return Constant(self.const * other.const)

        return super().__mul__(other)

    def __int__(self):
        """Return int(self)."""
        return int(self.const)

    def __float__(self):
        """Return float(self)."""
        return float(self.const)

    def __complex__(self):
        """Return complex(self)."""
        return complex(self.const)

    def __repr__(self):
        """Return repr(self)."""
        return "Constant({0!r})".format(self.const)
