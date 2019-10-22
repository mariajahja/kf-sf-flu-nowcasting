## Kalman Filter and Sensor Fusion for Flu Nowcasting
Code used to produce influenza nowcasts in the following paper. 

Jahja, M., Farrow, David C., Rosenfeld, R., Tibshirani, R.J. 
*Kalman Filter, Sensor Fusion, and Constrained Regression: Equivalences 
and Insights.* To appear: Neural Information Processing Systems (NeurIPS), 2019.

## Commands
The simulation in the paper was run by the following commands:
```sh
# neurips_main.py <starting epiweek> <ending epiweek> <name for output file>
> python3 neurips_main.py 201345 201420 1314
> python3 neurips_main.py 201445 201520 1415
> python3 neurips_main.py 201545 201620 1516
> python3 neurips_main.py 201645 201720 1617
> python3 neurips_main.py 201745 201820 1718
```
Note we do not produce predictions for the off-season. 

Trouble-shooting: Please ensure that the `src` module can be found on `PYTHONPATH`. A simple
workaround is to add the follow lines to the top of the simulation script:
```python
import sys
sys.path.append("~/<path>/kf-sf-flu-nowcasting")
```

## Considerations
The full nowcasting system uses several sources and sensors which are not 
made available to the public. Moreover, for historical reasons, several state 
signals (the influenza-like illness (ILI) response) are also kept private. 
For these reasons, the full simulation is not able to be publicly reproduced 
(please contact the first author for inquiries). 

A simplified set-up is given in `config.py`, which can be run to illustrate 
almost all the methods (except random forest trained on raw sources, which is
difficult to construct using only public information. Except for the change in 
inputs, this method follows exactly the same methodology as random forest 
trained on sensors). An example run is:
```sh
# neurips_main.py <starting epiweek> <ending epiweek> <name for output file>
> python3 neurips_main.py 201745 201820 1718
```
where, since the sensor features are sparse and small, there may be numerical
instability at very small values of regularization.


## Acknowledgements
Many of the utilities and frameworks were sourced and/or modified from CMU DELPHI:
https://github.com/cmu-delphi/. 
