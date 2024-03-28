import time
from selenium import webdriver
import subprocess
import json
import threading
import os
import sys


# Kill the server if a dash app is already running.
def kill_server():
    subprocess.run("lsof -t -i tcp:8080 | xargs kill -9", shell=True)  # run shell=True if on Linux


# Start Dash app.
def start_dash_app_frozen():
    path_dir = str(os.path.dirname(sys.executable))
    subprocess.Popen('dashboard', shell=False)


# Function to save browser session info.
def save_browser_session(input_driver):
    executor_url = input_driver.command_executor._url
    session_id = input_driver.session_id
    with open('browser_session.txt', 'w') as session_file:
        session_file.write(f'{str(executor_url)}\n{str(session_id)}')


# Start the web driver.
def start_driver():
    driver = webdriver.Chrome()
    time.sleep(5)
    driver.get('https://0.0.0.0:8080')
    save_browser_session(driver)


# Infinite while loop to keep the server running
def keep_server_running():
    while True:
        time.sleep(60)


def main():
    #kill_server()
    threading.Thread(target=start_dash_app_frozen).start()
    start_driver()
    keep_server_running()


if __name__ == '__main__':
    main()
