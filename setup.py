from setuptools import setup

setup(
    name="CellBro",
    version="0.1",
    long_description="",
    packages=['app'],
    include_package_data=True,
    zip_safe=False,
    install_requires=["Flask", "libsass >= 0.6.0", "Flask-JSGlue", "Jinja2"],
    setup_requires=['libsass >= 0.6.0'],
    sass_manifests={
        'app': ('static/sass', 'static/css', '/static/css')
    }
)