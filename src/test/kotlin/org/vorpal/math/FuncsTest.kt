// FuncsTest.kt
// By Sebastian Raaphorst, 2025.

package org.vorpal.math

import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.vorpal.org.vorpal.math.extendedGcd
import org.vorpal.org.vorpal.math.lcm
import org.vorpal.org.vorpal.math.powMod
import java.math.BigInteger
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals

class FuncsTest {
    @Test
    fun `Int 5^15 (mod 500)`() {
        assertEquals(125, powMod(5, 15, 500))
    }

    @Test
    fun `Long 10^20 (mod 600)`() {
        assertEquals(400L, powMod(10L, 20L, 600L))
    }

    @Test
    fun `BigInteger 99^1370 (mod 1234)`() {
        val x = BigInteger.valueOf(99)
        val y = BigInteger.valueOf(1370)
        val n = BigInteger.valueOf(1234)
        assertEquals(BigInteger.valueOf(445), powMod(x, y, n))
    }

    @Test
    fun `BigInteger 3^218 (mod 1000)`() {
        val x = BigInteger.valueOf(3)
        val y = BigInteger.valueOf(218)
        val n = BigInteger.valueOf(1000)
        assertEquals(BigInteger.valueOf(489), powMod(x, y, n))
        assertNotEquals(BigInteger.valueOf(1), powMod(x, y, n))
    }

    @Test
    fun `Int 1^1 (mod 0) throws ArithmeticException`() {
        assertThrows<ArithmeticException> { powMod(1, 1, 0) }
    }

    @Test
    fun `Long gcd and lcm of 32727312 and 3980289498`() {
        val a = 32_727_312L
        val b = 3_980_289_498L
        val (g, x, y) = extendedGcd(a, b)
        assertEquals(174L, g)
        assertEquals(174L, x * a + y * b)
        assertEquals(748_644_691_099_824L, lcm(a, b))
    }
}