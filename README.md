# THDSI
Total Harmonic Distortion for Speech Intelligibility

This is code I generated while completeing my PhD. Similar to Speech Transmission Index (STI), it requires both the STFT of a clean (original) signal and the STFT of a noisy (resulting) signal. It then computes the Total Harmonic distortion vs time array. Typically, Total harmonic distortion (THD) comuted on a transient signals by Head Acoustics Artemis or similar assumes the frequency bin with the highest amplitude to be the fundamental, but that is not always the case. THDSI uses the clean signal to determine the fundamental location and then computes THD on noisy signal. See chapter 5 my dissertation for a more detailed explanation: https://digitalcommons.mtu.edu/
