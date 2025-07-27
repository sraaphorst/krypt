// Fields.kt
// By Sebastian Raaphorst, 2025.

package org.vorpal.org.vorpal.math

import java.math.BigInteger

interface CommutativeRing<E> {
    val zero: E
    val one: E

    operator fun E.plus(b: E): E
    operator fun E.unaryMinus(): E
    operator fun E.minus(b: E): E = this + (-b)
    operator fun E.times(b: E): E
    fun characteristic(): BigInteger

    // Ways to call from within the ring: needed for implementations.
    fun add(a: E, b: E) = a + b
    fun subtract(a: E, b: E) = a - b
    fun multiply(a: E, b: E) = a * b
    fun negate(a: E) = -a
}

interface Field<E>: CommutativeRing<E> {
    fun E.inv(): E
    operator fun E.div(b: E): E = this * b.inv()
    fun invert(e: E): E = e.inv()
}

class PrimeField(val p: BigInteger) : Field<PrimeField.Element> {
    // For smaller prime finite fields.
    constructor(p: Int)  : this(BigInteger.valueOf(p.toLong()))
    constructor(p: Long) : this(BigInteger.valueOf(p))

    init {
        assert(p.isProbablePrime(100))
    }

    inner class Element(val n: BigInteger) {
        constructor(n: Int) : this(BigInteger.valueOf(n.toLong()))
        constructor(n: Long) : this(BigInteger.valueOf(n))

        init {
            require(n >= BigInteger.ZERO)
            require(n < p)
        }

        override fun toString() = n.toString()

        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (other !is Element) return false
            return n == other.n
        }

        override fun hashCode() = n.hashCode()
    }

    override val zero = Element(BigInteger.ZERO)
    override val one = Element(BigInteger.ONE)

    // To get an element of the field.
    operator fun invoke(n: BigInteger) = Element(n.mod(p))
    operator fun invoke(n: Int) = invoke(BigInteger.valueOf(n.toLong()))
    operator fun invoke(n: Long) = invoke(BigInteger.valueOf(n))

    override operator fun Element.plus(b: Element) = Element((n + b.n).mod(p))
    override operator fun Element.unaryMinus() = Element(n.negate().mod(p))
    override operator fun Element.times(b: Element) = Element((n * b.n).mod(p))
    override fun Element.inv(): Element =
        Element(powMod(n, p - BigInteger.TWO, p))

    override fun characteristic() = p

    fun generator(): Element {
        val phi = p - BigInteger.ONE

        val factors: Set<BigInteger> = run {
            tailrec fun factor(m: BigInteger, f: BigInteger, acc: Set<BigInteger>): Set<BigInteger> {
                return when {
                    f * f > m -> if (m > BigInteger.ONE) acc + m else acc
                    m % f == BigInteger.ZERO -> factor(m / f, f, acc + f)
                    else -> factor(m, f + BigInteger.ONE, acc)
                }
            }
            factor(phi, BigInteger.TWO, emptySet())
        }

        // Find the smallest candidate c >= 2 such that for every prime divisor q of phi:
        // c^(phi/q) mod p != 1.
        tailrec fun findRoot(c: BigInteger): BigInteger {
            if (c >= p) error("no generator found for F_$p")
            val isRoot = factors.all { q -> powMod(c, phi / q, p) != BigInteger.ONE }
            return if (isRoot) c else findRoot(c + BigInteger.ONE)
        }

        return Element(findRoot(BigInteger.TWO))
    }
}

class PolynomialRing<FE>(internal val ring: CommutativeRing<FE>) : CommutativeRing<PolynomialRing<FE>.Poly> {
    override val zero = Poly(listOf(ring.zero))
    override val one  = Poly(listOf(ring.one))

    inner class Poly(val coeffs: List<FE>) {
        operator fun plus(b: Poly): Poly {
            val max = maxOf(coeffs.size, b.coeffs.size)
            return Poly(List(max) { i ->
                val c = coeffs.getOrNull(i) ?: ring.zero
                val d = b.coeffs.getOrNull(i) ?: ring.zero
                ring.add(c, d)
            })
        }

        operator fun minus(b: Poly): Poly {
            val max = maxOf(coeffs.size, b.coeffs.size)
            return Poly(List(max) { i ->
                val c = coeffs.getOrNull(i) ?: ring.zero
                val d = b.coeffs.getOrNull(i) ?: ring.zero
                ring.subtract(c, d)
            })
        }

        operator fun times(b: Poly): Poly {
            val result = MutableList(coeffs.size + b.coeffs.size - 1) { ring.zero }
            for ((i, ai) in coeffs.withIndex()) {
                for ((j, bj) in b.coeffs.withIndex()) {
                    result[i + j] = ring.add(result[i + j], ring.multiply(ai, bj))
                }
            }
            return Poly(result)
        }

        operator fun PolynomialRing<FE>.Poly.unaryMinus(): PolynomialRing<FE>.Poly =
            Poly(coeffs.map { ring.negate(it) })

        override fun toString(): String {
            // build up “a·x^i” terms, skipping zero‐coefficients
            val terms = coeffs
                .mapIndexed { i, coef ->
                    println("coef is $coef class ${coef!!::class.simpleName}, ring.zero is ${ring.zero} class ${ring.zero!!::class.simpleName}, eq: ${coef == ring.zero}")
                    when {
                        coef == ring.zero -> { println("None"); null }
                        i == 0            -> coef.toString()
                        i == 1            -> "${coef}x"
                        else              -> "${coef}x^$i"
                    }
                }
                .filterNotNull()

            return if (terms.isEmpty()) "0"
            else terms.joinToString(" + ")
        }
    }

    override operator fun Poly.plus(b: Poly): Poly = this + b
    override operator fun Poly.unaryMinus(): Poly = -this
    override operator fun Poly.minus(b: Poly): Poly = this - b
    override operator fun Poly.times(b: Poly): Poly = this * b

    override fun characteristic(): BigInteger =
        ring.characteristic()

    /** handy constructor: ring(1,2,3) ⇒ 1 + 2x + 3x² */
    operator fun invoke(vararg cs: FE) = Poly(cs.toList())

    /** the “x” generator: 0 + 1 * x */
    val x = Poly(listOf(ring.zero, ring.one))
}

/**
 * Field extension made by taking the polynomial ring over a field mod a monic irreducible polynomial.
 */
class FieldExtension<FE>(private val ring: PolynomialRing<FE>,
                         private val modPoly: PolynomialRing<FE>.Poly) : Field<FieldExtension<FE>.Elt> {
    private val baseRing: Field<FE> = ring.ring as Field<FE>

    inner class Elt(val p: PolynomialRing<FE>.Poly) {
        internal fun reduce(q: PolynomialRing<FE>.Poly): PolynomialRing<FE>.Poly {
            val r = q.coeffs.toMutableList()
            val mdeg = modPoly.coeffs.lastIndex
            while (r.size - 1 >= mdeg) {
                val coeff = r.last()
                val shift = r.lastIndex - mdeg
                for ((i, mc) in modPoly.coeffs.withIndex()) {
                    r[shift + i] = baseRing.subtract(r[shift + i], baseRing.multiply(coeff, mc))
                }

            }
            while (r.isNotEmpty() && r.last() == baseRing.zero) r.removeAt(r.lastIndex)
            return ring.Poly(r)
        }

        fun inv(): Elt {
            // Use extended Euclid in the polynomial ring to find u,v with up + vm = 1
            // and return u mod m.
            tailrec fun ext(
                r0: PolynomialRing<FE>.Poly, r1: PolynomialRing<FE>.Poly,
                s0: PolynomialRing<FE>.Poly, s1: PolynomialRing<FE>.Poly): PolynomialRing<FE>.Poly =
                if (r1.coeffs.all { it == baseRing.zero }) s0
                else {
                    // polynomial division q = r0 / r1, rem = r0 % r1
                    val (q, rem) = polyDivMod(r0, r1)
                    ext(r1, rem, s1, ring.subtract(s0, ring.multiply(q, s1)))
                }
            val invPoly = ext(p, modPoly, ring(baseRing.one), ring(baseRing.zero))
            return Elt(reduce(invPoly))
        }

        override fun toString(): String {
            // build up “a·x^i” terms, skipping zero‐coefficients
            val terms = p.coeffs
                .mapIndexed { i, coef ->
                    when {
                        coef == baseRing.zero -> null
                        i == 0                -> coef.toString()
                        i == 1                -> "${coef}x"
                        else                  -> "${coef}x^$i"
                    }
                }
                .filterNotNull()

            return if (terms.isEmpty()) "0"
            else terms.joinToString(" + ")
        }
    }

    override val zero = Elt(ring.zero)
    override val one = Elt(ring.one)

    override operator fun Elt.plus(b: Elt) = Elt(reduce(ring.add(p, b.p)))
    override operator fun Elt.times(b: Elt) = Elt(reduce(ring.multiply(p, b.p)))
    override operator fun Elt.unaryMinus() = Elt(ring.negate(p))

    override fun Elt.inv(): Elt = this.inv()

    /**
     * Divide `dividend` by `divisor` in ring.F[x], returning (quotient, remainder).
     * Requires that the leading coefficient of divisor is invertible in ring.F.
     */
    private fun polyDivMod(
        dividend: PolynomialRing<FE>.Poly,
        divisor:  PolynomialRing<FE>.Poly
    ): Pair<PolynomialRing<FE>.Poly, PolynomialRing<FE>.Poly> {
        // Extract the coefficient lists: u will end up the remainder.
        val u = dividend.coeffs.toMutableList()
        val v = divisor.coeffs

        require(v.isNotEmpty() && v.last() != baseRing.zero) {
            "polyDivMod: division by zero polynomial"
        }

        val m = u.lastIndex       // degree of dividend
        val n = v.lastIndex       // degree of divisor
        if (m < n) {
            // degree(dividend) < degree(divisor) ⇒ quotient = 0, remainder = dividend
            return ring.Poly(emptyList()) to dividend
        }

        // prepare quotient coefficients (degree = m-n)
        val q = MutableList(m - n + 1) { baseRing.zero }

        // precompute the inverse of the leading coefficient of v
        val invLead = baseRing.invert(v[n])

        // Long division loop, from the highest down to n.
        for (k in m downTo n) {
            // compute factor = u[k] / v[n]
            val factor = u[k].let { baseRing.multiply(it, invLead) }
            q[k - n] = factor

            // subtract factor * x^(k-n) * v from u
            for (j in 0..n) {
                val idx = k - n + j
                u[idx] = baseRing.subtract(u[idx], baseRing.multiply(factor, v[j]))
            }
        }

        // Strip the leading zeros from the remainder.
        val rCoeffs = u.take(n).toMutableList()
        while (rCoeffs.isNotEmpty() && rCoeffs.last() == baseRing.zero) {
            rCoeffs.removeAt(rCoeffs.lastIndex)
        }

        return ring.Poly(q) to ring.Poly(rCoeffs)
    }

    override fun characteristic(): BigInteger =
        baseRing.characteristic()
}

fun main() {
    val f103 = PrimeField(103)
    val f103x = PolynomialRing(f103)
    val x = f103x.x
    val xsq = with(f103x) {
        x * x

    }
    val x2 = x * x * x + xsq + xsq + f103x.one
    println(x2.toString())

}