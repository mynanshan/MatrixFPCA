## Variance decomposition (in-mask, registered systole)

| alignment   | normalised   |   mask_pixels |   between |   within |   mean_dynamics |   subject_dynamics |
|:------------|:-------------|--------------:|----------:|---------:|----------------:|-------------------:|
| none        | False        |          7571 |     0.788 |    0.212 |           0.013 |              0.2   |
| none        | True         |          7571 |     0.775 |    0.225 |           0.013 |              0.212 |
| rigid       | False        |          5295 |     0.752 |    0.248 |           0.03  |              0.218 |
| rigid       | True         |          5295 |     0.732 |    0.268 |           0.032 |              0.236 |
| similarity  | False        |          4185 |     0.716 |    0.284 |           0.044 |              0.24  |
| similarity  | True         |          4185 |     0.688 |    0.312 |           0.047 |              0.265 |

## EF prediction from PCA baselines (TEST)

| alignment   | features           |   k |    R2 |   MAE |   AUC_EF40 |   lambda_ |
|:------------|:-------------------|----:|------:|------:|-----------:|----------:|
| none        | ED frame (static)  |  10 | 0.009 | 9.042 |      0.607 |   316.228 |
| none        | ED frame (static)  |  20 | 0.013 | 9.01  |      0.627 |   316.228 |
| none        | ED frame (static)  |  50 | 0.038 | 8.89  |      0.672 |  1000     |
| none        | Time average       |  10 | 0.037 | 8.948 |      0.642 |   316.228 |
| none        | Time average       |  20 | 0.065 | 8.835 |      0.684 |     0.01  |
| none        | Time average       |  50 | 0.089 | 8.724 |      0.718 |   100     |
| none        | ES - ED difference |  10 | 0.223 | 8.185 |      0.838 |   100     |
| none        | ES - ED difference |  20 | 0.265 | 7.975 |      0.86  |   100     |
| none        | ES - ED difference |  50 | 0.299 | 7.793 |      0.863 |   316.228 |
| none        | Stacked video      |  10 | 0.039 | 8.933 |      0.645 |   316.228 |
| none        | Stacked video      |  20 | 0.071 | 8.783 |      0.686 |     0.01  |
| none        | Stacked video      |  50 | 0.185 | 8.262 |      0.795 |   316.228 |
| rigid       | ED frame (static)  |  10 | 0.059 | 8.917 |      0.702 |     0.01  |
| rigid       | ED frame (static)  |  20 | 0.084 | 8.878 |      0.72  |     0.01  |
| rigid       | ED frame (static)  |  50 | 0.078 | 8.866 |      0.721 |   316.228 |
| rigid       | Time average       |  10 | 0.094 | 8.846 |      0.728 |    31.623 |
| rigid       | Time average       |  20 | 0.163 | 8.537 |      0.763 |     0.01  |
| rigid       | Time average       |  50 | 0.161 | 8.563 |      0.78  |   100     |
| rigid       | ES - ED difference |  10 | 0.26  | 8.002 |      0.872 |   100     |
| rigid       | ES - ED difference |  20 | 0.319 | 7.671 |      0.895 |   100     |
| rigid       | ES - ED difference |  50 | 0.349 | 7.686 |      0.91  |   100     |
| rigid       | Stacked video      |  10 | 0.099 | 8.813 |      0.736 |    31.623 |
| rigid       | Stacked video      |  20 | 0.187 | 8.47  |      0.783 |    31.623 |
| rigid       | Stacked video      |  50 | 0.384 | 7.383 |      0.904 |     0.01  |
| similarity  | ED frame (static)  |  10 | 0.061 | 8.886 |      0.694 |     0.01  |
| similarity  | ED frame (static)  |  20 | 0.107 | 8.772 |      0.73  |   100     |
| similarity  | ED frame (static)  |  50 | 0.119 | 8.752 |      0.746 |   100     |
| similarity  | Time average       |  10 | 0.149 | 8.627 |      0.771 |   100     |
| similarity  | Time average       |  20 | 0.202 | 8.41  |      0.798 |   316.228 |
| similarity  | Time average       |  50 | 0.27  | 8.064 |      0.829 |   100     |
| similarity  | ES - ED difference |  10 | 0.327 | 7.731 |      0.896 |   100     |
| similarity  | ES - ED difference |  20 | 0.35  | 7.611 |      0.913 |   100     |
| similarity  | ES - ED difference |  50 | 0.402 | 7.339 |      0.93  |   100     |
| similarity  | Stacked video      |  10 | 0.165 | 8.557 |      0.782 |   100     |
| similarity  | Stacked video      |  20 | 0.263 | 8.227 |      0.846 |   100     |
| similarity  | Stacked video      |  50 | 0.466 | 6.946 |      0.917 |    10     |
