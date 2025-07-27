// Funcs.kt
// By Sebastian Raaphorst, 2025.

package org.vorpal.org.vorpal.math

import java.math.BigInteger

/**
 * We want some general operators for all integer-like types to be able to calculate the gcd.
 * Since Kotlin does not have a Number interface, we create a type class to be able to run
 * the gcd and related algorithms.
 */
internal sealed interface EuclidOps<T> {
    fun zero(): T
    fun one(): T
    fun abs(x: T): T
    fun rem(a: T, b: T): T
    fun div(a: T, b: T): T
    fun sub(a: T, b: T): T
    fun mul(a: T, b: T): T
    fun eq(a: T, b: T): Boolean
    fun rmb(a: T): Boolean
    fun shr(a: T): T
}

internal object IntOps : EuclidOps<Int> {
    override fun zero(): Int = 0
    override fun one(): Int = 1
    override fun abs(x: Int): Int = kotlin.math.abs(x)
    override fun rem(a: Int, b: Int): Int = a % b
    override fun div(a: Int, b: Int): Int = a / b
    override fun sub(a: Int, b: Int): Int = a - b
    override fun mul(a: Int, b: Int): Int = a * b
    override fun eq(a: Int, b: Int): Boolean = a == b
    override fun rmb(a: Int): Boolean = (a and 1) == 1
    override fun shr(a: Int): Int = a shr 1
}

internal object LongOps : EuclidOps<Long> {
    override fun zero(): Long = 0L
    override fun one(): Long = 1
    override fun abs(x: Long): Long = kotlin.math.abs(x)
    override fun rem(a: Long, b: Long): Long = a % b
    override fun div(a: Long, b: Long): Long = a / b
    override fun sub(a: Long, b: Long): Long = a - b
    override fun mul(a: Long, b: Long): Long = a * b
    override fun eq(a: Long, b: Long): Boolean = a == b
    override fun rmb(a: Long): Boolean = (a and 1L) == 1L
    override fun shr(a: Long): Long = a shr 1
}

internal object BigIntegerOps : EuclidOps<BigInteger> {
    override fun zero(): BigInteger = BigInteger.ZERO
    override fun one(): BigInteger = BigInteger.ONE
    override fun abs(x: BigInteger): BigInteger = x.abs()
    override fun rem(a: BigInteger, b: BigInteger): BigInteger = a % b
    override fun div(a: BigInteger, b: BigInteger): BigInteger = a / b
    override fun sub(a: BigInteger, b: BigInteger): BigInteger = a - b
    override fun mul(a: BigInteger, b: BigInteger): BigInteger = a * b
    override fun eq(a: BigInteger, b: BigInteger): Boolean = a == b
    override fun rmb(a: BigInteger): Boolean = a.testBit(0)
    override fun shr(a: BigInteger): BigInteger = a.shiftRight(1)
}

internal fun <T> gcd(a: T, b: T, ops: EuclidOps<T>): T {
    tailrec fun aux(x: T = ops.abs(a), y: T = ops.abs(b)): T {
        return if (ops.eq(y, ops.zero())) x
        else aux(y, ops.rem(x, y))
    }
    return aux()
}

fun gcd(a: Int, b: Int): Int = gcd(a, b, IntOps)
fun gcd(a: Long, b: Long): Long = gcd(a, b, LongOps)
fun gcd(a: BigInteger, b: BigInteger): BigInteger = gcd(a, b, BigIntegerOps)

// To avoid large multiplications, divide the larger number by the gcd before multiplying.
fun lcm(a: Int, b: Int): Int = maxOf(a, b) / gcd(a, b, IntOps) * minOf(a, b)
fun lcm(a: Long, b: Long): Long = maxOf(a, b) / gcd(a, b, LongOps) * minOf(a, b)
fun lcm(a: BigInteger, b: BigInteger): BigInteger = maxOf(a, b) / gcd(a, b, BigIntegerOps) * minOf(a, b)

// Values for extended GCD output.
data class GcdValues<out T>(val gcd: T, val x: T, val y: T)

internal fun <T> extendedGcd(a0: T, b0: T, ops: EuclidOps<T>): GcdValues<T> {
    tailrec fun aux(r0: T, r1: T, s0: T, s1: T, t0: T, t1: T): GcdValues<T> {
        return if (ops.eq(r1, ops.zero())) {
            // r0 is gcd; s0 * a0 + t0 * b0 = r0
            GcdValues(r0, s0, t0)
        } else {
            val q  = ops.div(r0, r1)
            val r2 = ops.rem(r0, r1)
            val s2 = ops.sub(s0, ops.mul(q, s1))
            val t2 = ops.sub(t0, ops.mul(q, t1))

            aux(r1,r2, s1, s2, t1, t2)
        }
    }

    // initial values:
    // r0 = |a|, r1 = |b|
    // s0 = 1, s1 = 0, so s0 * a + s1 * b = a
    // t0 = 0, t1 = 1, so t0 * a + t1 * b = b
    return aux(
        ops.abs(a0),ops.abs(b0),
        ops.one(),ops.zero(),
        ops.zero(),ops.one()
    )
}

fun extendedGcd(a: Int, b: Int): GcdValues<Int> = extendedGcd(a, b, IntOps)
fun extendedGcd(a: Long, b: Long): GcdValues<Long> = extendedGcd(a, b, LongOps)
fun extendedGcd(a: BigInteger, b: BigInteger): GcdValues<BigInteger> = extendedGcd(a, b, BigIntegerOps)

/**
 * Calculate x^y (mod n) using fast multiplication.
 */
internal fun <T> powMod(x: T, y: T, n: T, ops: EuclidOps<T>): T {
    if (n == ops.zero())
        throw ArithmeticException("modulo by zero attempted")
    if (x == ops.one() || y == ops.zero())
        return ops.one()

    // At call i starting from i = 0, curr represents x^{2^i}.
    // If the least significant bit of remaining is 1, then we multiply total by curr.
    // We then square curr, shift remaining right by one position, and continue.
    tailrec fun aux(remaining: T = y,
                    curr: T = x,
                    total: T = ops.one()): T {
        if (ops.eq(remaining, ops.zero())) return total
        val newTotal = if (ops.rmb(remaining)) ops.rem(ops.mul(total, curr), n) else total
        val newCurr = ops.mul(curr, curr)
        val newRemaining = ops.shr(remaining)
        return aux(newRemaining, newCurr, newTotal)
    }
    return aux()
}

fun powMod(a: Int, b: Int, n: Int): Int = powMod(a, b, n, IntOps)
fun powMod(a: Long, b: Long, n: Long): Long = powMod(a, b, n, LongOps)
fun powMod(a: BigInteger, b: BigInteger, n: BigInteger): BigInteger = powMod(a, b, n, BigIntegerOps)
