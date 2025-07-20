package org.vorpal.org.vorpal.math

import java.math.BigInteger

sealed class GaloisField {
}

class FP(p: BigInteger) : GaloisField() {
    init {
        assert(p > BigInteger.ONE && p.isProbablePrime(100))
    }
}

fun main() {
    (1..63).forEach {
        val b = BigInteger.valueOf((2L shl it) - 1)
        println("$it -> $b is probably prime: ${b.isProbablePrime(100)}")
    }
}