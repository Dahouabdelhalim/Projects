# Install required packages
install.packages("rjson")
install.packages("jsonlite")
install.packages("httr")
install.packages("writexl")

# Load the packages.
library("rjson")
library("jsonlite")
library("httr")
library("writexl")

# Make a vector of the timeseries precinct data files
all_jsons <- c("https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-10-31T22:14:05.096Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-03T23:40:53.085Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:11:04.659Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:14:43.743Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:23:31.696Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:29:22.735Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:32:55.655Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:36:15.190Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:40:38.045Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:44:14.837Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:49:02.518Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T00:56:34.218Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:11:35.765Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:15:08.661Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:20:20.753Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:23:55.615Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:27:12.356Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:30:29.565Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:34:00.886Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:38:52.336Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:45:42.861Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:49:04.134Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:52:34.456Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T01:56:18.161Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:01:58.313Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:05:29.854Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:08:55.234Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:12:37.171Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:19:26.324Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:29:42.004Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:33:11.219Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:41:14.449Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:46:06.870Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:49:39.620Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:52:55.599Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T02:56:32.639Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:03:24.419Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:06:54.312Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:13:51.629Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:17:19.148Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:20:34.505Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:24:10.344Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:27:42.166Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:30:59.136Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:34:33.240Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:44:50.959Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:48:14.202Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:51:55.686Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:55:18.934Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T03:58:34.653Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:02:03.338Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:12:27.373Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:15:39.425Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:20:37.325Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:27:23.343Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:34:58.304Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:38:29.842Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:41:59.353Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:50:23.440Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:53:38.018Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T04:57:08.375Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T05:24:05.704Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T05:30:57.333Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T05:37:29.716Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T05:57:49.547Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T06:22:35.850Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T06:32:43.820Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T06:36:11.798Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T07:07:24.747Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T07:54:09.490Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T08:00:51.135Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T16:04:20.496Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T16:34:36.899Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T16:38:01.441Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T17:15:13.003Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T17:49:16.581Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T18:33:42.439Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T18:37:00.257Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T19:14:29.968Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T19:33:11.813Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T19:36:30.548Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T19:43:33.346Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T20:13:43.265Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T20:33:48.294Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T20:40:26.820Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T21:00:55.364Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T21:36:30.623Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T22:03:20.144Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T22:24:19.487Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T22:57:07.964Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-04T23:50:14.678Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T00:16:45.271Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T00:49:19.907Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T01:12:58.564Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T01:34:29.088Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T02:00:33.912Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T02:18:47.641Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T03:07:55.519Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T04:24:13.870Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T05:21:51.596Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T06:07:41.776Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T10:52:28.118Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T14:01:18.864Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T14:49:55.350Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T16:49:53.431Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T18:32:08.054Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T19:00:31.633Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T19:25:26.705Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T19:40:41.812Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T20:28:11.145Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T20:43:46.995Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T21:14:29.887Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T22:07:25.227Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T22:19:37.037Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T22:38:24.906Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T23:03:25.333Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-05T23:54:45.478Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T00:22:37.169Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T01:55:44.581Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T02:33:12.517Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T03:19:39.122Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T04:09:36.360Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T06:30:07.181Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T07:38:26.297Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T08:38:03.578Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T09:18:42.027Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T10:06:29.648Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T13:06:57.354Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T14:20:32.293Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T15:34:35.526Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T16:24:44.007Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T16:30:48.498Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T19:19:09.095Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T19:40:45.897Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T20:23:01.484Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T20:38:34.335Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T22:07:17.544Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T22:22:34.568Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T22:38:06.794Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T22:53:31.036Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T23:42:15.171Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T23:45:17.767Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-06T23:57:39.240Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T01:28:54.872Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T02:53:44.150Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T05:38:17.006Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T07:46:31.838Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T18:15:34.503Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T19:20:27.474Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T19:46:17.236Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-07T23:53:18.153Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-08T07:23:21.206Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-08T13:09:39.235Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-08T19:38:53.915Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T14:35:24.523Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T16:02:35.204Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T16:26:44.754Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T17:16:59.669Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T19:02:32.336Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T19:37:20.471Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T21:29:34.470Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T22:12:11.656Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T22:46:43.786Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-09T23:35:12.237Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T01:52:22.197Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T14:53:53.873Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T16:59:42.737Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T18:45:51.675Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T20:10:26.030Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T21:18:12.376Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-10T23:28:22.595Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-11T00:12:28.188Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-11T19:09:28.636Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-11T21:10:30.230Z.json",
                   "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/precincts/GAGeneral-2020-11-11T22:32:34.439Z.json")

# Load each timeseries precinct data file from the list and turn each file into a dataframe with just the abesntee votes from Fulton County
count <- 0
for (url in all_jsons) {
  count = count+1
  resp <- httr::GET(url)
  text_json <- httr::content(resp, type = 'text', encoding = "UTF-8")
  jfile <- fromJSON(text_json)
  df <- subset(as.data.frame(jfile), precincts.locality_name == "Fulton" & precincts.vote_type == "absentee")
  assign(paste("timepoint_", count, sep=""), df$precincts.results)
}

# Make a new matrix that shows the cummulative totals for each candidate at the first timepoint
timepoint <- 1
Trump <- sum(timepoint_1$trumpd)
Biden <- sum(timepoint_1$bidenj)
Jorgensen <- sum(timepoint_1$jorgensenj)
totals <- rbind(c(timepoint, Trump, Biden, Jorgensen))

# Name the columns
colnames(totals)<- c("Timepoint", "Trump","Biden","Jorgensen")

# Now add a rows to the matrix with cummulative totals for each candidate at each timepoint
timepoint <- 2
for (val in 1:174) {
  Trump <- sum(eval(parse(text=paste(paste("timepoint_", timepoint, sep = ""), "$trumpd", sep = ""))))
  Biden <- sum(eval(parse(text=paste(paste("timepoint_", timepoint, sep = ""), "$bidenj", sep = ""))))
  Jorgensen <- sum(eval(parse(text=paste(paste("timepoint_", timepoint, sep = ""), "$jorgensenj", sep = ""))))
  totals <- rbind(totals, c(timepoint, Trump, Biden, Jorgensen))
  timepoint = timepoint+1
}

# Convert the matrix to a dataframe
totals <- as.data.frame(totals)

# Export the dataframe as a Microsoft Excel file
write_xlsx(totals,"Fulton Absentee Votes.xlsx")
