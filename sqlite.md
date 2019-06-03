# Sqlite short introduction

[TOC]

###### create a database

use **sqlite3** in terminal to create a table

```sqlite
$ sqlite3 filename.db
```

###### create a table

use standard SQL language 

```sqlite
sqlite> create table tablename(column1 int, column2 text);
```

| Column1 | Column2 |
| ------- | ------- |
|         |         |

list tables

```sqlite
sqlite> .tables
```

view table details 

```sqlite
sqlite> .schema tablename
```

###### insert values into a table

```sqlite
sqlite> insert into tablename values(1, 'str1');
sqlite> insert into tablename values(2, 'str2');
```

| Column1 | Column2 |
| ------- | ------- |
| 1       | Str1    |
| 2       | Str2    |

###### select  data in a table

```sqlite
sqlite> select * from tablename
```

display headers

```sqlite
sqlite> .headers on
```

show settings

```sqlite
sqlite> .show
```

use **where**

```sqlite
sqlite> select * from tablename where column2 = 'str2';
sqlite> select * from tablename where column1 < 1;
```

select a column

```sqlite
sqlite> select column1 from tablename;
```

use **order by**

```sqlite
sqlite> select * from tablename where column1=1  order by column2;
```

use **limit** set how many output records 

```sqlite
sqlite> select * from tablename where column1=1 limit 3;
```

###### add a column

```sqlite
sqlite> alter table tablename add column column3 text;
```

###### rename a table

```sqlite
sqlite> alter table tablename rename to newname
```

###### update a table

update all rows if no **where** constraints

```sqlite
sqlite> update tablename set column1 = newvalue where constraints
```

###### delete data

delete all rows if no **where** constraints

```sqlite
sqlite> delete from tablename where constraints
```

###### drop a table

```sqlite
sqlite> drop table tablename
```

###### detach a database

```sqlite
sqlite> deteach database filename.db
```

###### import from csv

```sqlite
sqlite> .mode csv
sqlite> .import /path/csvname.csv csvtablename
```

###### Sqlite3 python

```python
import sqlite3
# connect to database
connection = sqlite3.connect('database.db')
cursor = connection.cursor()
# execute a sqlcommand
sqlcommand = 'SELECT col1,col2 FROM tablename WHERE constraints'
cursor.execute(sqlcommand)
# extract data as numpy array
data = cursor.fetchall()
# close
connection.close()
```

