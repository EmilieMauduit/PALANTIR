import pytest
import numpy as np

from palantir.prediction_tools import Star

@pytest.fixture
def star_data():
    """Fixture for creating a sample star with default values."""
    radius_data = {"models": ["Tout"], "radius": np.nan}
    return {
        "name": "TestStar",
        "mass": 1.0,
        "radius": radius_data,
        "age": 4.6,
        "obs_dist" : 1.0,
        "sp_type": "GV",
    }


@pytest.fixture
def star_instance(star_data):
    """Create a default Star instance for testing."""
    return Star(**star_data)


# Test for the Constructor and Property Getters
def test_constructor(star_instance):
    assert star_instance.name == "TestStar"
    assert star_instance.main_id == "TestStar"
    assert star_instance.mass == 1.0
    assert star_instance.age == 4.6e9
    assert star_instance.obs_dist == 1.0
    assert star_instance.sp_type == "GV"
    assert isinstance(star_instance.radius, float)
    assert isinstance(star_instance.luminosity, float)
    assert isinstance(star_instance.sp_type_code, int)
    assert isinstance(star_instance.rotperiod, float)



def test_effective_temperature(star_instance):
    assert star_instance.effective_temperature is None
    star_instance.compute_effective_temperature(value=10000.)
    assert star_instance.effective_temperature == 10000.
    star_instance.effective_temperature = None
    star_instance.compute_effective_temperature(value=np.nan)
    assert star_instance.effective_temperature > 0
    assert isinstance(star_instance.effective_temperature,float)

def test_magnetic_field(star_instance):
    star_instance.compute_effective_temperature(value=np.nan)
    
    value = {'model' : ['Bstar_original'], 'mag_field' : np.nan}
    star_instance.compute_magnetic_field(value)
    assert star_instance.magfield > 0
    assert isinstance(star_instance.magfield, float)


    value = {'model' : ['Bstar_polyfit'], 'mag_field' : np.nan}
    star_instance.magfield = None
    star_instance.compute_magnetic_field(value)
    assert star_instance.magfield > 0
    assert isinstance(star_instance.magfield, float)

    value = {'model' : ['Bstar_polyfit'], 'mag_field' : 12.}
    star_instance.magfield = None
    star_instance.compute_magnetic_field(value)
    assert star_instance.magfield == pytest.approx(12e-4, rel=1e-5)
    assert isinstance(star_instance.magfield, float)

    with pytest.raises(ValueError):
        star_instance.magfield = None
        star_instance.compute_magnetic_field({'model' : ['Bstar_original','Bstar_polyfit'], 'mag_field' : np.nan})