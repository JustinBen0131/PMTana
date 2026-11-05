# Calibrating the sPHENIX electromagnetic calorimeter with LED scans

A change in a calorimeter channel’s bias should produce a predictable change in its signal. Testing that prediction requires a stable light pulse. A single LED pulse width can leave some interface boards near noise and drive others close to saturation. Justin Bennett developed the C++/ROOT analysis here to choose a useful operating point for each board and measure how its response changes with bias. The sPHENIX EMCal channel map covers 24,576 towers connected to 384 interface boards.

This was a joint gain-adjustment study with Apurva Narde. Bennett’s code characterizes the response from the scans. The team then tested whether applying new bias settings produced the intended gain changes. Narde’s separate collaboration macro implements that comparison using the measured new-to-default response ratio.

## From a scan to a gain check

1. Associate towers with their sector and board, and associate each LED run with its pulse width and bias offset.
2. Rank measured pulse widths by how close the average board response is to 3,000 ADC counts at zero offset. Fit the response versus offset to an exponential and record its amplitude, slope and fit uncertainties. The number of boards with usable measurements depends on the scan.
3. For the joint check, compare the observed ratio after the adjustment with the intended gain, accounting for hardware limits. The gain targets were inputs to this work.

## Code and evidence

- [Board-level pulse selection and response fitting](PMTanalysisTestPulse.cpp). Start with `buildBestPulseWidthMap`. `produceCompare32vsBest` compares the chosen widths with one fixed width.
- [April 11 run settings](04-11-25-run-info.csv), [April 26 run settings](04-26-25-run-info.csv) and [tower-to-board map](vop.csv). Re-running the fit also requires the corresponding ROOT measurement files.
- [Pulse widths recorded by the collaboration](https://github.com/sPHENIX-Collaboration/analysis/blob/91c1ef7e63a48a3d277ffcc4882cd5d1c5bb9fc3/EMCal-Bias-Study/files/04-26-25-pulse-width-info-slim.csv) and [Narde’s measured-versus-intended response calculation](https://github.com/sPHENIX-Collaboration/analysis/blob/edea3062579cd18cd17b9ec489aa34b8174992a2/EMCal-Bias-Study/macros/analyze-new-response.C#L243-L273).
