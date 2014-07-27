path = 'C:\Users\Daniel Tavares\Desktop\daniel\sysident_data_2014.07.17\raw';

npts_exp = 10000;
npts_interval = 1000;
offset = -0.8/60/24;

start_date_array = {...
    '2014/07/17 19:46:00', ...
    '2014/07/17 19:48:30', ...
    '2014/07/17 19:51:01', ...
    '2014/07/17 19:53:32', ...
    '2014/07/17 19:56:04', ...
    '2014/07/17 19:58:34', ...
    '2014/07/17 20:01:04', ...
    '2014/07/17 20:03:35', ...
    '2014/07/17 20:06:06', ...
    '2014/07/17 20:08:38', ...
    '2014/07/17 20:11:08', ...
    '2014/07/17 20:13:38', ...
    '2014/07/17 20:16:09', ...
    '2014/07/17 20:18:40', ...
    '2014/07/17 20:21:10', ...
    '2014/07/17 20:34:06', ...
    '2014/07/17 20:36:41', ...
    '2014/07/17 20:39:12', ...
    '2014/07/17 20:41:43', ...
    '2014/07/17 20:44:15', ...
    '2014/07/17 20:46:45', ...
    '2014/07/17 20:49:15', ...
    '2014/07/17 20:51:45', ...
    '2014/07/17 20:54:16', ...
    '2014/07/17 20:56:48', ...
    '2014/07/17 20:59:18', ...
    '2014/07/17 21:01:48', ...
    '2014/07/17 21:04:18', ...
    '2014/07/17 21:06:49', ...
    '2014/07/17 21:09:20', ...
    '2014/07/17 21:25:48', ...
    '2014/07/17 21:28:18', ...
    '2014/07/17 21:30:48', ...
    '2014/07/17 21:33:19', ...
    '2014/07/17 21:33:54', ...
    '2014/07/17 21:36:25', ...
    '2014/07/17 21:38:55', ...
    '2014/07/17 21:41:26', ...
    '2014/07/17 21:43:57', ...
    '2014/07/17 21:46:29', ...
    '2014/07/17 21:49:00', ...
    '2014/07/17 21:51:30', ...
    '2014/07/17 21:54:01', ...
    '2014/07/17 21:56:32', ...
    '2014/07/17 21:59:04', ...
    '2014/07/17 22:09:28', ...
    '2014/07/17 22:11:58', ...
    '2014/07/17 22:14:29', ...
    '2014/07/17 22:17:03', ...
    '2014/07/17 22:19:37', ...
    };

idsaveexps(path, start_date_array, npts_exp, npts_interval, offset);