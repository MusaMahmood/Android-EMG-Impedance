package com.yeolabgt.mahmoodms.emgimpedance

import com.google.common.primitives.Doubles

/**
 * Created by mmahmood31 on 9/19/2017.
 * For Handling BLE incoming data packets.
 */

internal class DataChannel(var chEnabled: Boolean, MSBFirst: Boolean, //Classification:
                           private var classificationBufferSize: Int) {
    var dataBufferDoubles: DoubleArray? = null
    var packetCounter: Short = 0
    var totalDataPointsReceived: Int = 0
    var classificationBuffer: DoubleArray
    private var classificationBufferFloats: FloatArray

    init {
        this.packetCounter = 0
        this.totalDataPointsReceived = 0
        this.classificationBuffer = DoubleArray(classificationBufferSize)
        this.classificationBufferFloats = FloatArray(classificationBufferSize)
        Companion.MSBFirst = MSBFirst
    }

    /**
     * If 'dataBufferDoubles' is not null, concatenate new data using Guava lib
     * else: initialize dataBufferDoubles with new data.
     *
     * @param newDataPacket new data packet received via BLE>
     */
    fun handleNewData(newDataPacket: ByteArray) {
        this.totalDataPointsReceived += newDataPacket.size / 3
        val tempDoubleArray = DoubleArray(newDataPacket.size/3)
        for (i in 0 until newDataPacket.size / 3) {
            tempDoubleArray[i] = bytesToDouble(newDataPacket[3 * i], newDataPacket[3 * i + 1], newDataPacket[3 * i + 2])
        }
        addToDoublesBuffer(tempDoubleArray)
        if (this.dataBufferDoubles != null) {
            this.dataBufferDoubles = Doubles.concat(this.dataBufferDoubles, tempDoubleArray)
        } else {
            this.dataBufferDoubles = tempDoubleArray
        }
        this.packetCounter++
    }

    private fun addToDoublesBuffer(dataBufferDouble: DoubleArray) {
        if (this.classificationBufferSize > 0) {
            val newDataPoints = dataBufferDouble.size
            // Shift Data Backwards by N Amount
            System.arraycopy(this.classificationBuffer, newDataPoints, this.classificationBuffer, 0, this.classificationBufferSize - newDataPoints)
            // Copy new data to front of data
            dataBufferDouble.copyInto(this.classificationBuffer, this.classificationBufferSize - newDataPoints, 0, newDataPoints)
        }
    }

    fun resetBuffer() {
        this.dataBufferDoubles = null
        this.packetCounter = 0
    }

    companion object {
        private var MSBFirst: Boolean = false

        fun bytesToDoubleMPUAccel(a1: Byte, a2: Byte): Double {
            val unsigned: Int = unsignedBytesToInt(a1, a2, MSBFirst)
            return unsignedToSigned16bit(unsigned).toDouble() / 32767.0 * 16.0
        }

        fun bytesToDoubleMPUGyro(a1: Byte, a2: Byte): Double {
            val unsigned: Int = unsignedBytesToInt(a1, a2, MSBFirst)
            return unsignedToSigned16bit(unsigned).toDouble() / 32767.0 * 4000.0
        }

        fun bytesToFloat32(a1: Byte, a2: Byte, a3: Byte): Float {
            val unsigned = unsignedBytesToInt(a1, a2, a3, MSBFirst)
            return unsignedToSigned24bit(unsigned).toFloat() / 8388607.0.toFloat() * 2.25.toFloat()
        }

        fun bytesToFloat32(a1: Byte, a2: Byte): Float {
            val unsigned = unsignedBytesToInt(a1, a2, MSBFirst)
            return unsignedToSigned16bit(unsigned).toFloat() / 32767.0.toFloat() * 2.25.toFloat()
        }

        fun bytesToDouble(a1: Byte, a2: Byte): Double {
            val unsigned = unsignedBytesToInt(a1, a2, MSBFirst)
            return unsignedToSigned16bit(unsigned).toDouble() / 32767.0 * 2.25
        }


        fun bytesToDouble(a1: Byte, a2: Byte, a3: Byte): Double {
            val unsigned = unsignedBytesToInt(a1, a2, a3, MSBFirst)
            return unsignedToSigned24bit(unsigned).toDouble() / 8388607.0 * 2.25
        }

        private fun unsignedToSigned16bit(unsigned: Int): Int {
            return if (unsigned and 0x8000 != 0)
                -1 * (0x8000 - (unsigned and 0x8000 - 1))
            else
                unsigned
        }

        private fun unsignedToSigned24bit(unsigned: Int): Int {
            return if (unsigned and 0x800000 != 0) -1 * (0x800000 - (unsigned and 0x800000 - 1))
            else unsigned
        }

        private fun unsignedBytesToInt(b0: Byte, b1: Byte, MSBFirst: Boolean): Int {
            return if (MSBFirst)
                (unsignedByteToInt(b0) shl 8) + unsignedByteToInt(b1)
            else
                unsignedByteToInt(b0) + (unsignedByteToInt(b1) shl 8)
        }

        private fun unsignedBytesToInt(b0: Byte, b1: Byte, b2: Byte, MSBFirst: Boolean): Int {
            return if (MSBFirst)
                (unsignedByteToInt(b0) shl 16) + (unsignedByteToInt(b1) shl 8) + unsignedByteToInt(b2)
            else
                unsignedByteToInt(b0) + (unsignedByteToInt(b1) shl 8) + (unsignedByteToInt(b2) shl 16)
        }

        private fun unsignedByteToInt(b: Byte): Int {
            return (b.toInt() and 0xFF)
        }

//        private fun unsignedToSigned(unsignedInt: Int, size: Int): Int {
//            var unsigned = unsignedInt
//            if (unsigned and (1 shl size - 1) != 0) unsigned = -1 * ((1 shl size - 1) - (unsigned and (1 shl size - 1) - 1))
//            return unsigned
//        }
    }
}
