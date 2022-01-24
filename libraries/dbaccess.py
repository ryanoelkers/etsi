""" This class is used for basic functions not spcecific to this code base such as:
 logging, file writing, and testing."""
from config import Configuration
from libraries.utils import Utils
import psycopg2
import pandas as pd
import os


class DBaccess:

    @staticmethod
    def get_ticid_query(path, ticids):
        """ This function will query the tic for all columns based on the TICID.

        :parameter path - The location of the SQL file
        :parameter ticids - The ticid of the stars

        :return A string with the SQL
        """

        # read in teh appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        tics = ', '.join(ticids)
        sql_cmd = sql_cmd.replace('%(ticid)s', tics)

        return sql_cmd

    @staticmethod
    def get_starlist_query(path, cen_ra, cen_de, mx_dist):
        """ This function will replace the determined centers and distance in the SQL FILE and then return a string
        useful for querying the data base.

        :parameter path - The location of the SQL file
        :parameter cen_ra - The center RA coordinate
        :parameter cen_de - The center Dec coordinate
        :parameter mx_dist - The maximum edge to center distance

        :return A string with the SQL
        """
        # read in the appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        sql_cmd = sql_cmd.replace('%(cen_ra)s', str(cen_ra))
        sql_cmd = sql_cmd.replace('%(cen_de)s', str(cen_de))
        sql_cmd = sql_cmd.replace('%(mx_dist)s', str(mx_dist))

        return sql_cmd

    @staticmethod
    def query_tic7_bulk(sql_cmd, stassunlab):
        """ This function will query tic v 7 on tessdev, and return a data frame based on the query used. This function
        will confirm a dumped .csv file does not exist prior to accessing the database.

        :parameter sql_cmd - The sql query to use
        :parameter stassunlab - The current computer being used for the program

        :return df - a data frame with the query results
        """
        # set up the connection object based on whether you are on tessdev
        if stassunlab != 'tessdev':
            conn = psycopg2.connect(host="129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        if stassunlab == 'tessdev':
            conn = psycopg2.connect(host="129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        # set up the cursor object
        cur = conn.cursor()

        Utils.log("Querying TICv7 for the next " + str(Configuration.BULK_QUERY) + " stars.",
                  "info", Configuration.LOG_SCREEN)

        # generate the data frame with the queried results
        df = pd.read_sql_query(sql_cmd, conn)

        Utils.log("Query complete!", "info", Configuration.LOG_SCREEN)

        # shut it down
        cur.close()
        conn.close()

        return df

    @staticmethod
    def query_tic7(sql_cmd, out_path, file_name, stassunlab):
        """ This function will query tic v 7 on tessdev, and return a data frame based on the query used. This function
        will confirm a dumped .csv file does not exist prior to accessing the database.

        :parameter sql_cmd - The sql query to use
        :parameter out_path - The output for the file
        :parameter file_name - The desired filename
        :parameter stassunlab - The current computer being used for the program

        :return df - a data frame with the query results
        """
        # set up the connection object based on whether you are on tessdev
        if stassunlab != 'tessdev':
            conn = psycopg2.connect(host="129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        if stassunlab == 'tessdev':
            conn = psycopg2.connect(host="129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        # set up the cursor object
        cur = conn.cursor()

        if os.path.isfile(out_path + file_name) == 1:
            Utils.log("Legacy file found, not querying TICv7.", "info")
            # read in from a file
            df = pd.read_csv(out_path + file_name, index_col=0)
            Utils.log("CSV read complete.", "info")

        if os.path.isfile(out_path + file_name) == 0:
            Utils.log("Querying TICv7...", "info")
            # generate the data frame with the queried results
            df = pd.read_sql_query(sql_cmd, conn)

            Utils.log("Query complete, dumping to .csv file " + out_path + file_name, "info")
            # dump the file to csv
            df.to_csv(out_path + file_name)

            Utils.log("Dump complete.", "info")

        # convert the pk name to a TICID
        df = df.rename(columns={'pk': 'TICID'})

        # shut it down
        cur.close()
        conn.close()

        return df

