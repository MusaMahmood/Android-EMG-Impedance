//
// Created by mahmoodms on 4/3/2017.
//

#include "rt_nonfinite.h"
#include "ssvep_filter_f32.h"
#include "downsample_250Hz.h"
#include "ecg_bandstop_250Hz.h"
#include "extract_fftz.h"

/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define  LOG_TAG "jniExecutor-cpp"
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

extern "C" {
JNIEXPORT jdoubleArray JNICALL
Java_com_yeolabgt_mahmoodms_emgimpedance_DeviceControlActivity_jextractFFT(JNIEnv *env,
                                                                           jobject jobject1,
                                                                           jdoubleArray data) {
    jdouble *X = env->GetDoubleArrayElements(data, nullptr);
    if (X == nullptr) LOGE("ERROR - X is null");
    double Y[800];
    jdoubleArray m_result = env->NewDoubleArray(800);
    extract_fftz(X, Y);
    env->SetDoubleArrayRegion(m_result, 0, 800, Y);
    return m_result;
}
}

extern "C" {
JNIEXPORT jdoubleArray JNICALL
Java_com_yeolabgt_mahmoodms_emgimpedance_DeviceControlActivity_jecgBandStopFilter(
        JNIEnv *env, jobject jobject1, jdoubleArray data) {
    jdouble *X1 = env->GetDoubleArrayElements(data, nullptr);
    double Y[1000]; // First two values = Y; last 499 = cPSD
    if (X1 == nullptr) LOGE("ERROR - C_ARRAY IS nullptr");
    jdoubleArray m_result = env->NewDoubleArray(1000);
    ecg_bandstop_250Hz(X1, Y);
    env->SetDoubleArrayRegion(m_result, 0, 1000, Y);
    return m_result;
}
}


extern "C" {
JNIEXPORT jdoubleArray JNICALL
Java_com_yeolabgt_mahmoodms_emgimpedance_DeviceControlActivity_jdownSample(
        JNIEnv *env, jobject jobject1, jdoubleArray data, jint Fs) {
    jdouble *X1 = env->GetDoubleArrayElements(data, nullptr);
    int Xsize[1] = {Fs * 4};
    double Y[1000]; // First two values = Y; last 499 = cPSD
    int Ysize[2]; // First two values = Y; last 499 = cPSD
    if (X1 == nullptr) LOGE("ERROR - C_ARRAY IS nullptr");
    jdoubleArray m_result = env->NewDoubleArray(1000);
    downsample_250Hz(X1, Xsize, Fs, &Y[0], Ysize);
    env->SetDoubleArrayRegion(m_result, 0, 1000, Y);
    return m_result;
}
}


extern "C" {
JNIEXPORT jfloatArray JNICALL
Java_com_yeolabgt_mahmoodms_emgimpedance_DeviceControlActivity_jSSVEPCfilter(
        JNIEnv *env, jobject jobject1, jdoubleArray data) {
    jdouble *X1 = env->GetDoubleArrayElements(data, nullptr);
    float Y[1000]; // First two values = Y; last 499 = cPSD
    if (X1 == nullptr) LOGE("ERROR - C_ARRAY IS nullptr");
    jfloatArray m_result = env->NewFloatArray(1000);
    ssvep_filter_f32(X1, Y);
    env->SetFloatArrayRegion(m_result, 0, 1000, Y);
    return m_result;
}
}

extern "C" {
JNIEXPORT jdoubleArray JNICALL
/**
 *
 * @param env
 * @param jobject1
 * @return array of frequencies (Hz) corresponding to a raw input signal.
 */
Java_com_yeolabgt_mahmoodms_emgimpedance_DeviceControlActivity_jLoadfFFT(
        JNIEnv *env, jobject jobject1) {
    jdoubleArray m_result = env->NewDoubleArray(800);
    double fFFT[800];
    for (int i = 0; i < 800; i++) {
        fFFT[i] = (double) 900.0 + (i * 0.2500);
    }
    env->SetDoubleArrayRegion(m_result, 0, 800, fFFT);
    return m_result;
}
}

extern "C" {
JNIEXPORT jint JNICALL
Java_com_yeolabgt_mahmoodms_emgimpedance_DeviceControlActivity_jmainInitialization(
        JNIEnv *env, jobject obj, jboolean initialize) {
    if (!(bool) initialize) {
        downsample_250Hz_initialize();
        ecg_bandstop_250Hz_initialize();
        extract_fftz_initialize();
        return 0;
    } else {
        return -1;
    }
}
}
