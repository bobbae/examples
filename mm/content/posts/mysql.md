---
title:  "Installing mysql on Ubuntu 18.04 wsl2"
date: 2020-05-20

#header:
#   teaser: "/assets/images/500x300.png" 
tags: 
  - "mysql"
---

After installing with the usual

```
$ sudo apt install mysql-server
$ sudo service mysql start
$ sudo mysql_secure_installation
```

I *could not* login by doing

```
$ mysql -u root
ERROR 1698 (28000): Access denied for user 'root'@'localhost'
```


Looking at the authentication plugin shows:

```
mysql> SELECT user,authentication_string,plugin,host FROM mysql.user;
+------------------+-------------------------------------------+-----------------------+----
-------+
| user             | authentication_string                     | plugin                | hos
t      |
+------------------+-------------------------------------------+-----------------------+----
-------+
| root             | *2470C0C06DEE42FD1618BB99005ADCA2EC9D1E19 | auth_socket           | loc
alhost |
| mysql.session    | *THISISNOTAVALIDPASSWORDTHATCANBEUSEDHERE | mysql_native_password | loc
alhost |
| mysql.sys        | *THISISNOTAVALIDPASSWORDTHATCANBEUSEDHERE | mysql_native_password | loc
alhost |
| debian-sys-maint | *0F2B3822681E14406843203F5708B1CD54AC57C8 | mysql_native_password | loc
alhost |
+------------------+-------------------------------------------+-----------------------+----
-------+
4 rows in set (0.00 sec)
```


So I tried

```
$ sudo mysql
mysql> ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'password';
ERROR 1819 (HY000): Your password does not satisfy the current policy requirements
```

which does not work until you first uninstall validate_password plugin for some reason.

```
mysql> uninstall plugin validate_password;
Query OK, 0 rows affected (0.02 sec)
```

After validate_password is removed, you can

```
mysql> ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'password';
mysql> FLUSH PRIVILEGES;
Query OK, 0 rows affected (0.00 sec)

mysql> \q
Bye
```

Now you can do

```
$ mysql -u root -p
Enter password:
Welcome to the MySQL monitor.  Commands end with ; or \g.
Your MySQL connection id is 7
Server version: 5.7.30-0ubuntu0.18.04.1 (Ubuntu)

Copyright (c) 2000, 2020, Oracle and/or its affiliates. All rights reserved.

Oracle is a registered trademark of Oracle Corporation and/or its
affiliates. Other names may be trademarks of their respective
owners.

Type 'help;' or '\h' for help. Type '\c' to clear the current input statement.

mysql>
```


