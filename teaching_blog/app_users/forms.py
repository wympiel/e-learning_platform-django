from django import forms
from django.contrib.auth.models import User
from .models import UserProfileInfo
from django.contrib.auth.forms import UserCreationForm


class UserForm(UserCreationForm):
    email = forms.EmailField()

    class Meta():
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2')

        labels = {
            'password1': 'Hasło',
            'password2': 'Potwierdź hasło'
        }


class UserProfileInfoForm(forms.ModelForm):
    bio = forms.CharField(required=False)
    teacher = 'trener'
    student = 'student'
    parent = 'rodzic'
    user_types = [
        (student, 'student'),
        (parent, 'rodzic')
    ]
    user_type = forms.ChoiceField(required=True, choices=user_types)

    class Meta():
        model = UserProfileInfo
        fields = ('bio', 'profile_pic', 'user_type')
